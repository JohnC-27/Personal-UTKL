#!/usr/bin/env python3
"""
optimize_kde.py — Fit an adaptive ROOT.TKDE template to a target TH1 by minimizing chi-squared.

Optimizes global bandwidth scale (rho) and amplitude scale (alpha) using scipy.optimize.minimize.
Uses the unweighted TKDE constructor; kernel and iteration mode are set via the option string.
"""

from __future__ import annotations

import array
import math
import os
import sys
from dataclasses import dataclass
from typing import Callable, Sequence, Tuple

import ROOT
from scipy.optimize import minimize

ROOT.gErrorIgnoreLevel = ROOT.kWarning

# -----------------------------------------------------------------------------
# Input / output configuration (edit for your workflow)
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
OPT_METHOD = "Nelder-Mead"  # or "L-BFGS-B" with bounds below
INITIAL_RHO = 0.2
INITIAL_ALPHA = 1.0
LBFGS_BOUNDS = ((0.05, None), (1e-9, None))  # rho > 0.05, alpha > 0


@dataclass
class KdeFitContext:
  """Immutable inputs shared by the loss function and optimizers."""

  data: array.array
  x_min: float
  x_max: float
  target: ROOT.TH1
  options: str
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


def histogram_to_unweighted_sample(hist: ROOT.TH1) -> array.array:
  """
  Expand non-empty TH1 bins into an unweighted event sample for TKDE.

  Each bin with content > 0 contributes int(round(content)) copies of its bin center.
  Bins with fractional content still contribute at least one event when content >= 0.5.
  """
  points: list[float] = []
  for i in range(1, hist.GetNbinsX() + 1):
    content = hist.GetBinContent(i)
    if content <= 0:
      continue
    n_rep = max(1, int(round(content)))
    center = hist.GetBinCenter(i)
    points.extend([center] * n_rep)

  if not points:
    raise ValueError("histogram has no positive bin content; cannot build TKDE sample")

  return array.array("d", points)


def make_tkde(rho: float, ctx: KdeFitContext) -> ROOT.TKDE:
  """Construct TKDE with the required unweighted C++ binding signature."""
  return ROOT.TKDE(
    len(ctx.data),
    ctx.data,
    ctx.x_min,
    ctx.x_max,
    ctx.options,
    float(rho),
  )


def chi_squared_vs_hist(
  kde_func: ROOT.TF1,
  hist: ROOT.TH1,
  alpha: float,
) -> float:
  """
  Standard chi-squared sum over histogram bins.

  Skips bins with zero content or zero error to avoid division-by-zero.
  """
  chi2 = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    observed = hist.GetBinContent(i)
    if observed <= 0:
      continue
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    expected = alpha * kde_func.Eval(hist.GetBinCenter(i))
    diff = observed - expected
    chi2 += (diff * diff) / (err * err)
  return chi2


def build_chi_squared_loss(ctx: KdeFitContext) -> Callable[[Sequence[float]], float]:
  """
  Return a closure loss(params) for scipy, where params = [rho, alpha].

  Re-instantiates TKDE on every evaluation because rho changes the adaptive fit.
  """

  _last_kde: list[ROOT.TKDE] = []  # keep ROOT object alive across Eval calls

  def chi_squared_loss(params: Sequence[float]) -> float:
    rho, alpha = float(params[0]), float(params[1])
    if rho <= 0 or alpha <= 0:
      return math.inf

    kde = make_tkde(rho, ctx)
    kde_func = kde.GetFunction()
    if kde_func is None:
      return math.inf

    kde_func.SetNpx(ctx.npx)
    _last_kde.clear()
    _last_kde.append(kde)

    return chi_squared_vs_hist(kde_func, ctx.target, alpha)

  return chi_squared_loss


def count_chi2_degrees_of_freedom(hist: ROOT.TH1) -> int:
  """Number of bins included in the chi-squared sum (positive content and error)."""
  ndf = 0
  for i in range(1, hist.GetNbinsX() + 1):
    if hist.GetBinContent(i) > 0 and hist.GetBinError(i) > 0:
      ndf += 1
  return ndf


def optimize_kde_chi2(
  ctx: KdeFitContext,
  x0: Tuple[float, float] = (INITIAL_RHO, INITIAL_ALPHA),
  method: str = OPT_METHOD,
) -> Tuple[float, float, object]:
  """Run scipy.optimize.minimize and return (rho, alpha, result)."""
  loss = build_chi_squared_loss(ctx)
  x0_arr = [float(x0[0]), float(x0[1])]

  if method.upper() == "L-BFGS-B":
    result = minimize(
      loss,
      x0_arr,
      method="L-BFGS-B",
      bounds=LBFGS_BOUNDS,
      options={"maxiter": 500, "ftol": 1e-8},
    )
  else:
    result = minimize(
      loss,
      x0_arr,
      method="Nelder-Mead",
      options={"maxiter": 2000, "xatol": 1e-4, "fatol": 1e-4},
    )

  rho_opt, alpha_opt = float(result.x[0]), float(result.x[1])
  return rho_opt, alpha_opt, result


def scaled_kde_tf1(
  rho: float,
  alpha: float,
  ctx: KdeFitContext,
) -> Tuple[ROOT.TKDE, ROOT.TF1]:
  """Build final TKDE and a TF1 wrapper f(x) = alpha * kde(x)."""
  kde = make_tkde(rho, ctx)
  raw = kde.GetFunction()
  if raw is None:
    raise RuntimeError("TKDE::GetFunction() returned null")

  raw.SetNpx(ctx.npx)

  def scaled(x, _p):
    return alpha * raw.Eval(x[0])

  name = f"kde_fit_rho{rho:.4g}_a{alpha:.4g}"
  fit_fn = ROOT.TF1(name, scaled, ctx.x_min, ctx.x_max, 0)
  fit_fn.SetNpx(ctx.npx)
  fit_fn._hold_kde = kde
  fit_fn._hold_raw = raw
  return kde, fit_fn


def kde_template_histogram(
  kde_func: ROOT.TF1,
  ref_hist: ROOT.TH1,
  alpha: float,
  name: str = "kde_template",
) -> ROOT.TH1D:
  """Binned template: alpha * KDE evaluated at ref_hist bin centers (portable for I/O)."""
  nb = ref_hist.GetNbinsX()
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()
  out = ROOT.TH1D(name, "alpha #times KDE template;X;Expected", nb, xlo, xhi)
  out.SetDirectory(0)
  for i in range(1, nb + 1):
    out.SetBinContent(i, alpha * kde_func.Eval(ref_hist.GetBinCenter(i)))
  return out


def save_fit_results(
  outfile: str,
  target: ROOT.TH1,
  kde_shape: ROOT.TF1,
  template_hist: ROOT.TH1,
  rho: float,
  alpha: float,
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
  shape.SetDirectory(fout)
  shape.Write()

  template_hist.SetDirectory(fout)
  template_hist.Write()

  meta = ROOT.TNamed(
    "fit_meta",
    f"rho={rho};alpha={alpha};chi2={chi2};ndf={ndf};reduced_chi2={chi2 / max(ndf - 2, 1)}",
  )
  meta.Write()

  fout.Write()
  fout.Close()


def main() -> int:
  print(f"Reading {TARGET_HIST_NAME!r} from {INPUT_ROOT_FILE}")
  target = open_input_histogram(INPUT_ROOT_FILE, TARGET_HIST_NAME, PROJECTION_AXIS)

  data = histogram_to_unweighted_sample(target)
  x_min = target.GetXaxis().GetXmin()
  x_max = target.GetXaxis().GetXmax()

  ctx = KdeFitContext(
    data=data,
    x_min=x_min,
    x_max=x_max,
    target=target,
    options=TKDE_OPTIONS,
  )

  print(f"TKDE sample: {len(data)} unweighted points on [{x_min}, {x_max}]")
  print(f"Options: {TKDE_OPTIONS!r}")
  print(f"Optimizer: {OPT_METHOD}  initial (rho, alpha) = ({INITIAL_RHO}, {INITIAL_ALPHA})")

  rho_opt, alpha_opt, opt_result = optimize_kde_chi2(ctx)
  print(f"Optimization success={opt_result.success}  message: {opt_result.message}")
  print(f"Optimal rho = {rho_opt:.6g}  alpha = {alpha_opt:.6g}")

  kde, fit_fn = scaled_kde_tf1(rho_opt, alpha_opt, ctx)
  raw_shape = fit_fn._hold_raw
  chi2 = chi_squared_vs_hist(raw_shape, target, alpha_opt)
  ndf = count_chi2_degrees_of_freedom(target)
  npar = 2
  reduced = chi2 / max(ndf - npar, 1)

  print(f"chi2 = {chi2:.6g}  ndf = {ndf}  reduced chi2 = {reduced:.6g}")

  template = kde_template_histogram(raw_shape, target, alpha_opt)

  os.makedirs(os.path.dirname(os.path.abspath(OUTPUT_ROOT_FILE)), exist_ok=True)
  save_fit_results(
    OUTPUT_ROOT_FILE, target, raw_shape, template, rho_opt, alpha_opt, chi2, ndf
  )
  del kde, fit_fn  # release after Clone() in save_fit_results
  print(f"Wrote fit to {OUTPUT_ROOT_FILE}")
  return 0 if opt_result.success else 1


if __name__ == "__main__":
  sys.exit(main())
