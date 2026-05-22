#!/usr/bin/env python3
"""
optimize_kde.py — Fit a weighted adaptive ROOT.TKDE template to a target TH1 via chi-squared.

Minimizes chi2/ndf over bandwidth scale (rho) and x-axis scale (x_scale) with
Nelder-Mead; amplitude (alpha) is profiled analytically at each trial.
Model: alpha * kde(x_scale * x). Uses the weighted TKDE constructor.
"""

from __future__ import annotations

import array
import math
import os
import sys
from dataclasses import dataclass
from typing import Sequence, Tuple

import ROOT
from scipy.optimize import minimize

ROOT.gErrorIgnoreLevel = ROOT.kWarning

# -----------------------------------------------------------------------------
# Input / output configuration
# -----------------------------------------------------------------------------
INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "nominal.root"
)
TARGET_HIST_NAME = "nominalxyposMM1"  # TH2D; projected to X below
PROJECTION_AXIS = "x"  # "x" or "y" for TH2 projections; ignored for TH1

OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "kde_chi2_fit.root"
)

# TKDE option string (kernel type must live here, not as an integer argument)
TKDE_OPTIONS = "KernelType:Gaussian;Iteration:Adaptive;Mirror:noMirror"

# Optimizer defaults
OPT_METHOD = "Nelder-Mead"
NM_MAXITER = 2000
RHO_MIN = 0.11  # floor: blocks the near-zero spike / interpolating limit
INITIAL_RHO = 0.2
INITIAL_X_SCALE = 1.0
X_SCALE_MIN = 0.9
X_SCALE_MAX = 1.1
NPAR = 2  # rho, x_scale fitted; alpha profiled (linear WLS)


@dataclass
class KdeFitContext:
  """Immutable inputs shared by the loss function and optimizers."""

  data: array.array
  weights: array.array
  x_min: float
  x_max: float
  target: ROOT.TH1
  options: str
  ndf: int
  npx: int = 2000


def open_input_histogram(
  filepath: str,
  hist_name: str,
  projection_axis: str = "x",
) -> ROOT.TH1:
  """Open a ROOT file and return the target TH1 (or a 1D projection of a TH2)."""
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  obj = tfile.Get(hist_name)
  if not obj:
    tfile.Close()
    raise KeyError(f"object {hist_name!r} not found in {filepath}")

  if obj.InheritsFrom("TH2"):
    axis = projection_axis.lower()
    if axis == "x":
      hist = obj.ProjectionX(f"{hist_name}_px")
    elif axis == "y":
      hist = obj.ProjectionY(f"{hist_name}_py")
    else:
      tfile.Close()
      raise ValueError(f"projection_axis must be 'x' or 'y', got {projection_axis!r}")
  elif obj.InheritsFrom("TH1"):
    hist = obj
  else:
    tfile.Close()
    raise TypeError(f"{hist_name!r} is not TH1/TH2")

  hist.SetDirectory(0)
  tfile.Close()
  return hist


def histogram_to_weighted_sample(hist: ROOT.TH1) -> Tuple[array.array, array.array]:
  """
  Build weighted TKDE inputs from TH1 bin centers (data) and bin contents (weights).

  Only bins with content > 0 are included. data and weights have identical length.
  """
  xs: list[float] = []
  ws: list[float] = []
  for i in range(1, hist.GetNbinsX() + 1):
    content = hist.GetBinContent(i)
    if content <= 0:
      continue
    xs.append(hist.GetBinCenter(i))
    ws.append(content)

  if not xs:
    raise ValueError("histogram has no positive bin content; cannot build weighted TKDE")

  return array.array("d", xs), array.array("d", ws)


def make_tkde(rho: float, ctx: KdeFitContext) -> ROOT.TKDE:
  """Construct TKDE with the weighted C++ binding signature."""
  if len(ctx.data) != len(ctx.weights):
    raise ValueError("data and weights arrays must have identical length")
  return ROOT.TKDE(
    len(ctx.data),
    ctx.data,
    ctx.weights,
    ctx.x_min,
    ctx.x_max,
    ctx.options,
    float(rho),
  )


def scaled_kde_value(
  kde_func: ROOT.TF1,
  x: float,
  alpha: float,
  x_scale: float,
) -> float:
  """alpha * kde(x_scale * x)."""
  return alpha * kde_func.Eval(x_scale * x)


def optimal_alpha(
  kde_func: ROOT.TF1,
  hist: ROOT.TH1,
  x_scale: float,
) -> float:
  """Profiled amplitude minimizing chi2 at fixed kde shape and x_scale."""
  num, den = 0.0, 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    observed = hist.GetBinContent(i)
    shape = kde_func.Eval(x_scale * hist.GetBinCenter(i))
    w = 1.0 / (err * err)
    num += w * observed * shape
    den += w * shape * shape
  if den <= 0:
    return 1.0
  return num / den


def max_x_scale_for_domain(ctx: KdeFitContext) -> float:
  """Largest x_scale such that x_scale * bin_center stays in [x_min, x_max]."""
  limit = float("inf")
  for i in range(1, ctx.target.GetNbinsX() + 1):
    if ctx.target.GetBinError(i) <= 0:
      continue
    x = ctx.target.GetBinCenter(i)
    if x > 0:
      limit = min(limit, ctx.x_max / x)
    elif x < 0:
      limit = min(limit, ctx.x_min / x)
  return limit * 0.999 if math.isfinite(limit) else X_SCALE_MAX


def chi_squared_vs_hist(
  kde_func: ROOT.TF1,
  hist: ROOT.TH1,
  alpha: float,
  x_scale: float,
) -> float:
  """
  Standard chi-squared sum over histogram bins.

  Skips bins with zero error to avoid division-by-zero.
  """
  chi2 = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    observed = hist.GetBinContent(i)
    expected = scaled_kde_value(kde_func, hist.GetBinCenter(i), alpha, x_scale)
    diff = observed - expected
    chi2 += (diff * diff) / (err * err)
  return chi2


def count_chi2_bins(hist: ROOT.TH1) -> int:
  """Number of histogram bins entering the chi-squared sum (error > 0)."""
  n = 0
  for i in range(1, hist.GetNbinsX() + 1):
    if hist.GetBinError(i) > 0:
      n += 1
  return n


def degrees_of_freedom(nbins: int, npar: int = NPAR) -> int:
  """ndf = nbins - nparams."""
  return max(nbins - npar, 1)


def pack_optimizer_params(
  rho: float,
  x_scale: float,
  x_scale_max: float = X_SCALE_MAX,
) -> Tuple[float, float]:
  """Map physical (rho, x_scale) to unconstrained Nelder-Mead variables [u, w]."""
  rho = max(float(rho), RHO_MIN)
  x_hi = min(x_scale_max, X_SCALE_MAX)
  x_scale = min(max(float(x_scale), X_SCALE_MIN), x_hi)
  u = math.sqrt(rho - RHO_MIN)
  w = math.log(x_scale)
  return u, w


def unpack_optimizer_params(
  params: Sequence[float],
  x_scale_max: float = X_SCALE_MAX,
) -> Tuple[float, float]:
  """Map [u, w] to rho >= RHO_MIN and x_scale in [X_SCALE_MIN, x_scale_max]."""
  u, w = float(params[0]), float(params[1])
  x_hi = min(x_scale_max, X_SCALE_MAX)
  x_scale = min(max(math.exp(w), X_SCALE_MIN), x_hi)
  return RHO_MIN + u * u, x_scale


def evaluate_at_rho_xscale(
  rho: float,
  x_scale: float,
  ctx: KdeFitContext,
  hold: list[ROOT.TKDE] | None = None,
) -> Tuple[float, float, float]:
  """Return (chi2, reduced chi2, profiled alpha) for trial (rho, x_scale)."""
  kde = make_tkde(rho, ctx)
  kde_func = kde.GetFunction()
  if kde_func is None:
    return math.inf, math.inf, 1.0
  kde_func.SetNpx(ctx.npx)
  if hold is not None:
    hold.clear()
    hold.append(kde)
  alpha = optimal_alpha(kde_func, ctx.target, x_scale)
  chi2 = chi_squared_vs_hist(kde_func, ctx.target, alpha, x_scale)
  reduced = chi2 / max(ctx.ndf, 1)
  return chi2, reduced, alpha


def build_reduced_chi2_loss(ctx: KdeFitContext, x_scale_max: float):
  """
  Return loss(internal_params) = chi2/ndf for Nelder-Mead.

  internal_params = [u, w] with rho = RHO_MIN + u^2, x_scale = exp(w); alpha profiled.
  """
  _last_kde: list[ROOT.TKDE] = []

  def loss(params: Sequence[float]) -> float:
    rho, x_scale = unpack_optimizer_params(params, x_scale_max)
    _chi2, reduced, _alpha = evaluate_at_rho_xscale(rho, x_scale, ctx, _last_kde)
    return reduced

  return loss


def initial_fit_guess(ctx: KdeFitContext) -> Tuple[float, float, float]:
  """Starting (rho, profiled alpha, x_scale) at nominal bandwidth."""
  _chi2, _reduced, alpha = evaluate_at_rho_xscale(INITIAL_RHO, INITIAL_X_SCALE, ctx)
  return INITIAL_RHO, alpha, INITIAL_X_SCALE


def optimize_kde_chi2(
  ctx: KdeFitContext,
  x0: Tuple[float, float] | None = None,
) -> Tuple[float, float, float, object]:
  """Minimize chi2/ndf over rho and x_scale; alpha profiled at each evaluation."""
  x_scale_max = min(max_x_scale_for_domain(ctx), X_SCALE_MAX)
  if x0 is None:
    rho0, _alpha0, x_scale0 = initial_fit_guess(ctx)
    x0 = (rho0, x_scale0)
  loss = build_reduced_chi2_loss(ctx, x_scale_max)
  u0, w0 = pack_optimizer_params(x0[0], x0[1], x_scale_max)

  result = minimize(
    loss,
    [u0, w0],
    method="Nelder-Mead",
    options={"maxiter": NM_MAXITER, "xatol": 1e-4, "fatol": 1e-4},
  )

  rho_opt, x_scale_opt = unpack_optimizer_params(result.x, x_scale_max)
  _chi2, _reduced, alpha_opt = evaluate_at_rho_xscale(rho_opt, x_scale_opt, ctx)
  return rho_opt, alpha_opt, x_scale_opt, result


def scaled_kde_tf1(
  rho: float,
  alpha: float,
  x_scale: float,
  ctx: KdeFitContext,
) -> Tuple[ROOT.TKDE, ROOT.TF1]:
  """Build final weighted TKDE and a TF1 wrapper f(x) = alpha * kde(x_scale * x)."""
  kde = make_tkde(rho, ctx)
  raw = kde.GetFunction()
  if raw is None:
    raise RuntimeError("TKDE::GetFunction() returned null")

  raw.SetNpx(ctx.npx)

  def scaled(x, _p):
    return scaled_kde_value(raw, x[0], alpha, x_scale)

  name = f"kde_fit_rho{rho:.4g}_a{alpha:.4g}_s{x_scale:.4g}"
  fit_fn = ROOT.TF1(name, scaled, ctx.x_min, ctx.x_max, 0)
  fit_fn.SetNpx(ctx.npx)
  fit_fn._hold_kde = kde
  fit_fn._hold_raw = raw
  fit_fn._x_scale = x_scale
  return kde, fit_fn


def kde_template_histogram(
  kde_func: ROOT.TF1,
  ref_hist: ROOT.TH1,
  alpha: float,
  x_scale: float,
  name: str = "kde_template",
) -> ROOT.TH1D:
  """Binned template: alpha * KDE(x_scale * x) at ref_hist bin centers."""
  nb = ref_hist.GetNbinsX()
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()
  out = ROOT.TH1D(name, "#alpha#timesKDE(#it{s}#cdotx);X;Expected", nb, xlo, xhi)
  out.SetDirectory(0)
  for i in range(1, nb + 1):
    x = ref_hist.GetBinCenter(i)
    out.SetBinContent(i, scaled_kde_value(kde_func, x, alpha, x_scale))
  return out


def save_fit_results(
  outfile: str,
  target: ROOT.TH1,
  kde_shape: ROOT.TF1,
  template_hist: ROOT.TH1,
  rho: float,
  alpha: float,
  x_scale: float,
  chi2: float,
  ndf: int,
) -> None:
  """Write histograms, KDE shape TF1, binned template, and fit metadata."""
  fout = ROOT.TFile.Open(outfile, "RECREATE")
  if not fout or fout.IsZombie():
    raise OSError(f"cannot create output file: {outfile}")

  hist_out = target.Clone("target_hist")
  hist_out.SetDirectory(fout)
  hist_out.Write()

  shape = kde_shape.Clone("kde_shape")
  fout.cd()
  shape.Write()

  template_hist.SetDirectory(fout)
  template_hist.Write()

  reduced = chi2 / max(ndf, 1)
  meta = ROOT.TNamed(
    "fit_meta",
    f"rho={rho};alpha={alpha};x_scale={x_scale};chi2={chi2};ndf={ndf};reduced_chi2={reduced}",
  )
  fout.cd()
  meta.Write()

  fout.Write()
  fout.Close()


def main() -> int:
  print(f"Reading {TARGET_HIST_NAME!r} from {INPUT_ROOT_FILE}")
  target = open_input_histogram(INPUT_ROOT_FILE, TARGET_HIST_NAME, PROJECTION_AXIS)

  data, weights = histogram_to_weighted_sample(target)
  x_min = target.GetXaxis().GetXmin()
  x_max = target.GetXaxis().GetXmax()

  nbins = count_chi2_bins(target)
  ndf = degrees_of_freedom(nbins)

  ctx = KdeFitContext(
    data=data,
    weights=weights,
    x_min=x_min,
    x_max=x_max,
    target=target,
    options=TKDE_OPTIONS,
    ndf=ndf,
  )

  print(f"Weighted TKDE: {len(data)} points on [{x_min}, {x_max}]")
  print(f"Options: {TKDE_OPTIONS!r}")
  x_scale_max = min(max_x_scale_for_domain(ctx), X_SCALE_MAX)
  rho0, alpha0, x_scale0 = initial_fit_guess(ctx)
  print(
    f"rho >= {RHO_MIN}  x_scale in [{X_SCALE_MIN}, {x_scale_max:.4g}]  "
    f"ndf = {ndf} (nbins - {NPAR} fit params, alpha profiled)  "
    f"model: alpha*kde(x_scale*x)"
  )
  print(
    f"Optimizer: {OPT_METHOD} on (rho, x_scale); alpha profiled  "
    f"initial (rho, x_scale) = ({rho0:.4g}, {x_scale0:.4g})  "
    f"alpha@start = {alpha0:.4g}"
  )

  rho_opt, alpha_opt, x_scale_opt, opt_result = optimize_kde_chi2(ctx, (rho0, x_scale0))
  print(f"Optimization success={opt_result.success}  message: {opt_result.message}")
  print(
    f"Optimal rho = {rho_opt:.6g}  alpha = {alpha_opt:.6g}  x_scale = {x_scale_opt:.6g}"
  )

  kde, fit_fn = scaled_kde_tf1(rho_opt, alpha_opt, x_scale_opt, ctx)
  raw_shape = fit_fn._hold_raw
  chi2 = chi_squared_vs_hist(raw_shape, target, alpha_opt, x_scale_opt)
  reduced = chi2 / max(ndf, 1)

  print(f"chi2 = {chi2:.6g}  nbins = {nbins}  ndf = {ndf}  chi2/ndf = {reduced:.6g}")

  template = kde_template_histogram(raw_shape, target, alpha_opt, x_scale_opt)

  os.makedirs(os.path.dirname(os.path.abspath(OUTPUT_ROOT_FILE)), exist_ok=True)
  save_fit_results(
    OUTPUT_ROOT_FILE,
    target,
    raw_shape,
    template,
    rho_opt,
    alpha_opt,
    x_scale_opt,
    chi2,
    ndf,
  )
  del kde, fit_fn
  print(f"Wrote fit to {OUTPUT_ROOT_FILE}")
  return 0 if opt_result.success else 1


if __name__ == "__main__":
  sys.exit(main())
