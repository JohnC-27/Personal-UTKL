#!/usr/bin/env python3
"""
weighted_adaptive_kde.py — Weighted adaptive ROOT.TKDE bandwidth scan.

Scans bandwidth scale rho from RHO_SCAN_MIN to RHO_SCAN_MAX, profiles alpha at each
point, records chi-squared vs the target TH1, writes results to ROOT, and plots
chi2/ndf vs bandwidth.
"""
import argparse
import array
import math
import os
import sys
from dataclasses import dataclass
from typing import Iterator, List, Tuple

import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.SetBatch(True)

INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "nominal.root"
  # nominal.root or mz_nominal_2000bin_run1.root (or run2)
)
  # Capital N for Mo's hists. lower for nominal.root
TARGET_HIST_NAME = "nominalxyposMM1"  # TH2D; projected to X below
PROJECTION_AXIS = "y"  # "x" or "y" for TH2 projections; ignored for TH1

OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "weighted_fixed_kde.root"
)
CHI2_SCAN_PNG = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "weighted_fixed_kde_chi2_scan.png"
)

TKDE_BASE_OPTIONS = "KernelType:Gaussian;Iteration:Fixed"
TKDE_OPTIONS_NO_MIRROR = f"{TKDE_BASE_OPTIONS};Mirror:noMirror"
TKDE_OPTIONS_MIRROR = f"{TKDE_BASE_OPTIONS};Mirror:MirrorBoth"
TKDE_OPTIONS = TKDE_OPTIONS_MIRROR

# True: optimize w*unmirrored + (1-w)*mirrored. False: single KDE (TKDE_OPTIONS).
USE_LINEAR_COMBO = True

RHO_SCAN_MIN = 0.18
RHO_SCAN_MAX = 0.18
RHO_SCAN_STEP = 0.0001





@dataclass
class RhoScanPoint:
  rho: float
  chi2: float
  reduced_chi2: float
  alpha: float
  mix: float  # weight on unmirrored KDE; mirrored weight is 1 - mix


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
  npx: int = 10000


@dataclass
class DistributionStats:
  mean: float
  std: float
  integral: float


def rho_scan_values(
  rho_min: float = RHO_SCAN_MIN,
  rho_max: float = RHO_SCAN_MAX,
  step: float = RHO_SCAN_STEP,
) -> Iterator[float]:
  """Inclusive bandwidth grid from rho_min to rho_max."""
  rho = rho_min
  eps = 1e-12 * max(abs(rho_max), 1.0)
  while rho <= rho_max + eps:
    yield rho
    rho += step


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


def drop_edge_bins(hist: ROOT.TH1) -> ROOT.TH1:
  nb = hist.GetNbinsX()
  if nb < 3:
    raise ValueError("need at least 3 bins")

  xlo = hist.GetBinLowEdge(2)
  xhi = hist.GetBinLowEdge(nb)  # right edge of bin nb-1
  out = hist.Clone("hist_trimmed")
  out.SetDirectory(0)
  out.SetBins(nb - 2, xlo, xhi)

  for i in range(2, nb):
    out.SetBinContent(i - 1, hist.GetBinContent(i))
    out.SetBinError(i - 1, hist.GetBinError(i))

  if out.GetSumw2N() == 0 and hist.GetSumw2N() > 0:
    out.Sumw2()

  return out


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


def make_tkde(
  rho: float,
  ctx: KdeFitContext,
  *,
  options: str | None = None,
) -> ROOT.TKDE:
  """Construct TKDE with the weighted C++ binding signature."""
  if len(ctx.data) != len(ctx.weights):
    raise ValueError("data and weights arrays must have identical length")
  return ROOT.TKDE(
    len(ctx.data),
    ctx.data,
    ctx.weights,
    ctx.x_min,
    ctx.x_max,
    options if options is not None else ctx.options,
    float(rho),
  )


def combined_kde_shape(
  unmirrored: ROOT.TF1,
  mirrored: ROOT.TF1,
  mix: float,
  x: float,
) -> float:
  """Linear blend: mix * unmirrored + (1 - mix) * mirrored."""
  return mix * unmirrored.Eval(x) + (1.0 - mix) * mirrored.Eval(x)


def make_combined_kde_tf1(
  unmirrored: ROOT.TF1,
  mirrored: ROOT.TF1,
  mix: float,
  x_min: float,
  x_max: float,
  npx: int,
  name: str = "kde_combined_shape",
) -> ROOT.TF1:
  """TF1 wrapper for the unnormalized combined KDE shape."""

  def blended(x, _p):
    return combined_kde_shape(unmirrored, mirrored, mix, x[0])

  out = ROOT.TF1(name, blended, x_min, x_max, 0)
  out.SetNpx(npx)
  out._hold_unmirrored = unmirrored
  out._hold_mirrored = mirrored
  return out


def tkde_pair_at_rho(
  rho: float,
  ctx: KdeFitContext,
) -> Tuple[ROOT.TKDE, ROOT.TF1, ROOT.TKDE, ROOT.TF1]:
  """Build weighted fixed unmirrored and mirrored TKDEs at fixed rho."""
  kde_unmirrored = make_tkde(rho, ctx, options=TKDE_OPTIONS_NO_MIRROR)
  kde_mirrored = make_tkde(rho, ctx, options=TKDE_OPTIONS_MIRROR)

  unmirrored = kde_unmirrored.GetFunction()
  mirrored = kde_mirrored.GetFunction()
  if unmirrored is None or mirrored is None:
    raise RuntimeError("TKDE::GetFunction() returned null")

  unmirrored.SetNpx(ctx.npx)
  mirrored.SetNpx(ctx.npx)
  return kde_unmirrored, unmirrored, kde_mirrored, mirrored


def evaluate_at_mix(
  mix: float,
  unmirrored: ROOT.TF1,
  mirrored: ROOT.TF1,
  ctx: KdeFitContext,
) -> Tuple[float, float, float]:
  """Return (chi2, reduced chi2, profiled alpha) for a fixed mix weight."""
  shape = make_combined_kde_tf1(
    unmirrored,
    mirrored,
    mix,
    ctx.x_min,
    ctx.x_max,
    ctx.npx,
    name=f"kde_mix{mix:.4g}",
  )
  alpha = optimal_alpha(shape, ctx.target)
  chi2 = chi_squared_vs_hist(shape, ctx.target, alpha)
  reduced = chi2 / max(ctx.ndf, 1)
  return chi2, reduced, alpha


def _golden_section_minimize(
  func,
  lo: float,
  hi: float,
  *,
  tol: float = 1e-4,
  max_iter: int = 80,
) -> float:
  """Minimize a unimodal function on [lo, hi] via golden-section search."""
  phi = 0.5 * (math.sqrt(5.0) - 1.0)
  a, b = lo, hi
  c = b - phi * (b - a)
  d = a + phi * (b - a)
  fc = func(c)
  fd = func(d)

  for _ in range(max_iter):
    if b - a <= tol:
      break
    if fc < fd:
      b, d, fd = d, c, fc
      c = b - phi * (b - a)
      fc = func(c)
    else:
      a, c, fc = c, d, fd
      d = a + phi * (b - a)
      fd = func(d)

  return 0.5 * (a + b)


def optimize_mix_weight(
  unmirrored: ROOT.TF1,
  mirrored: ROOT.TF1,
  ctx: KdeFitContext,
) -> Tuple[float, float, float, float]:
  """
  Minimize chi2/ndf over mix in [0, 1] with alpha profiled at each trial.

  Returns (mix, alpha, chi2, reduced_chi2).
  """
  cache: dict[float, Tuple[float, float, float]] = {}

  def loss(mix: float) -> float:
    if mix not in cache:
      cache[mix] = evaluate_at_mix(mix, unmirrored, mirrored, ctx)
    _chi2, reduced, _alpha = cache[mix]
    if not math.isfinite(reduced):
      return math.inf
    return reduced

  mix_opt = _golden_section_minimize(loss, 0.0, 1.0)
  chi2, reduced, alpha = evaluate_at_mix(mix_opt, unmirrored, mirrored, ctx)
  return mix_opt, alpha, chi2, reduced


def optimal_alpha(kde_func: ROOT.TF1, hist: ROOT.TH1) -> float:
  """Profiled amplitude minimizing chi2 at fixed kde shape."""
  num, den = 0.0, 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    observed = hist.GetBinContent(i)
    shape = kde_func.Eval(hist.GetBinCenter(i))
    w = 1.0 / (err * err)
    num += w * observed * shape
    den += w * shape * shape
  if den <= 0:
    return 1.0
  return num / den


def chi_squared_vs_hist(
  kde_func: ROOT.TF1,
  hist: ROOT.TH1,
  alpha: float,
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
    expected = alpha * kde_func.Eval(hist.GetBinCenter(i))
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


def chi2_bin_contribution(
  hist: ROOT.TH1,
  kde_func: ROOT.TF1,
  alpha: float,
  bin_index: int,
) -> float | None:
  """Per-bin chi2 term; None if bin is skipped (zero error)."""
  err = hist.GetBinError(bin_index)
  if err <= 0:
    return None
  observed = hist.GetBinContent(bin_index)
  expected = alpha * kde_func.Eval(hist.GetBinCenter(bin_index))
  diff = observed - expected
  return (diff * diff) / (err * err)


def print_chi2_contributions(
  hist: ROOT.TH1,
  kde_func: ROOT.TF1,
  alpha: float,
  *,
  exclude_first_last: bool = True,
) -> None:
  """Print per-bin chi2 contributions and total for interior bins."""
  nb = hist.GetNbinsX()
  skip = {1, nb} if exclude_first_last else set()
  if exclude_first_last:
    scope = f"bins 2..{nb - 1} (excluding first and last)"
  else:
    scope = f"bins 1..{nb}"

  print(f"\nChi2 contributions ({scope}), alpha = {alpha:.6g}:")
  print(f"{'bin':>5} {'x [cm]':>10} {'data':>12} {'model':>12} {'chi2_i':>12}")

  total = 0.0
  n_used = 0
  for i in range(1, nb + 1):
    if i in skip:
      continue
    contrib = chi2_bin_contribution(hist, kde_func, alpha, i)
    if contrib is None:
      print(f"{i:5d}  (skipped: zero error)")
      continue
    observed = hist.GetBinContent(i)
    expected = alpha * kde_func.Eval(hist.GetBinCenter(i))
    total += contrib
    n_used += 1
    x = hist.GetBinCenter(i)
    print(
      f"{i:5d} {x:10.4f} {observed:12.4g} {expected:12.4g} {contrib:12.6g}"
    )

  ndf = degrees_of_freedom(n_used)
  reduced = total / ndf
  print(f"Sum chi2 = {total:.6g}  (n_bins = {n_used}, ndf = {ndf}, chi2/ndf = {reduced:.6g})")


def fit_nparams(use_linear_combo: bool = USE_LINEAR_COMBO) -> int:
  """Free parameters in the fit; alpha is always profiled."""
  return 2 if use_linear_combo else 1


def degrees_of_freedom(
  nbins: int,
  *,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> int:
  """ndf = nbins - nparams."""
  return max(nbins - fit_nparams(use_linear_combo), 1)


def evaluate_at_rho(
  rho: float,
  ctx: KdeFitContext,
  *,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> Tuple[float, float, float, float]:
  """
  Return (chi2, reduced chi2, profiled alpha, mix) at fixed rho.

  mix weights the unmirrored KDE when use_linear_combo is True; otherwise 0.
  """
  if not use_linear_combo:
    kde = make_tkde(rho, ctx)
    kde_func = kde.GetFunction()
    if kde_func is None:
      return math.inf, math.inf, 1.0, 0.0
    kde_func.SetNpx(ctx.npx)
    alpha = optimal_alpha(kde_func, ctx.target)
    chi2 = chi_squared_vs_hist(kde_func, ctx.target, alpha)
    reduced = chi2 / max(ctx.ndf, 1)
    return chi2, reduced, alpha, 0.0

  try:
    _kde_u, unmirrored, _kde_m, mirrored = tkde_pair_at_rho(rho, ctx)
  except RuntimeError:
    return math.inf, math.inf, 1.0, 0.0

  mix, alpha, chi2, reduced = optimize_mix_weight(unmirrored, mirrored, ctx)
  return chi2, reduced, alpha, mix


def scan_rho_bandwidths(
  ctx: KdeFitContext,
  *,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> List[RhoScanPoint]:
  """Evaluate chi2 at each rho on the scan grid."""
  results: List[RhoScanPoint] = []
  for rho in rho_scan_values():
    chi2, reduced, alpha, mix = evaluate_at_rho(
      rho,
      ctx,
      use_linear_combo=use_linear_combo,
    )
    if not math.isfinite(chi2):
      print(f"  rho={rho:.4g}: KDE failed, skipped")
      continue
    results.append(
      RhoScanPoint(
        rho=rho,
        chi2=chi2,
        reduced_chi2=reduced,
        alpha=alpha,
        mix=mix,
      )
    )
  if not results:
    raise RuntimeError("no successful KDE evaluations on bandwidth grid")
  return results


def scan_to_graphs(
  scan: List[RhoScanPoint],
) -> Tuple[ROOT.TGraph, ROOT.TGraph, ROOT.TGraph]:
  """Build TGraphs: chi2, chi2/ndf, and alpha vs rho."""
  n = len(scan)
  g_chi2 = ROOT.TGraph(n)
  g_chi2.SetName("chi2_vs_rho")
  g_chi2.SetTitle("#chi^{2} vs bandwidth #rho;#rho;#chi^{2}")

  g_reduced = ROOT.TGraph(n)
  g_reduced.SetName("reduced_chi2_vs_rho")
  g_reduced.SetTitle("#chi^{2}/ndf vs bandwidth #rho;#rho;#chi^{2}/ndf")

  g_alpha = ROOT.TGraph(n)
  g_alpha.SetName("alpha_vs_rho")
  g_alpha.SetTitle("#alpha vs bandwidth #rho;#rho;#alpha")

  for i, pt in enumerate(scan):
    g_chi2.SetPoint(i, pt.rho, pt.chi2)
    g_reduced.SetPoint(i, pt.rho, pt.reduced_chi2)
    g_alpha.SetPoint(i, pt.rho, pt.alpha)
  return g_chi2, g_reduced, g_alpha


def best_scan_point(scan: List[RhoScanPoint]) -> RhoScanPoint:
  """Point with minimum reduced chi-squared."""
  return min(scan, key=lambda p: p.reduced_chi2)


def scaled_kde_tf1(
  rho: float,
  alpha: float,
  ctx: KdeFitContext,
  *,
  mix: float = 1.0,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> Tuple[ROOT.TKDE | None, ROOT.TKDE | None, ROOT.TF1, ROOT.TF1]:
  """
  Build scaled KDE fit function at the best scan point.

  Returns (kde_unmirrored, kde_mirrored, shape, fit_fn). Single-KDE mode leaves
  kde_unmirrored None and uses kde_mirrored for the lone TKDE object.
  """
  if not use_linear_combo:
    kde = make_tkde(rho, ctx)
    raw = kde.GetFunction()
    if raw is None:
      raise RuntimeError("TKDE::GetFunction() returned null")

    raw.SetNpx(ctx.npx)

    def scaled(x, _p):
      return alpha * raw.Eval(x[0])

    name = f"kde_rho{rho:.4g}_a{alpha:.4g}"
    fit_fn = ROOT.TF1(name, scaled, ctx.x_min, ctx.x_max, 0)
    fit_fn.SetNpx(ctx.npx)
    fit_fn._hold_kde = kde
    fit_fn._hold_raw = raw
    return None, kde, raw, fit_fn

  kde_unmirrored, unmirrored, kde_mirrored, mirrored = tkde_pair_at_rho(rho, ctx)
  combined_shape = make_combined_kde_tf1(
    unmirrored,
    mirrored,
    mix,
    ctx.x_min,
    ctx.x_max,
    ctx.npx,
    name=f"kde_shape_rho{rho:.4g}_mix{mix:.4g}",
  )

  def scaled(x, _p):
    return alpha * combined_shape.Eval(x[0])

  name = f"kde_rho{rho:.4g}_mix{mix:.4g}_a{alpha:.4g}"
  fit_fn = ROOT.TF1(name, scaled, ctx.x_min, ctx.x_max, 0)
  fit_fn.SetNpx(ctx.npx)
  fit_fn._hold_kde_unmirrored = kde_unmirrored
  fit_fn._hold_kde_mirrored = kde_mirrored
  fit_fn._hold_combined_shape = combined_shape
  fit_fn._hold_unmirrored = unmirrored
  fit_fn._hold_mirrored = mirrored
  return kde_unmirrored, kde_mirrored, combined_shape, fit_fn


def kde_template_histogram(
  kde_func: ROOT.TF1,
  ref_hist: ROOT.TH1,
  alpha: float,
  name: str = "kde_template",
  x_axis_title: str = PROJECTION_AXIS
) -> ROOT.TH1D:
  """Binned template: alpha * KDE(x) at ref_hist bin centers."""
  nb = ref_hist.GetNbinsX()
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()

  if x_axis_title.lower() == 'x':
    x_axis_title = 'x [cm]'
  elif x_axis_title.lower() == 'y':
    x_axis_title = 'y [cm]'
  else:
    x_axis_title = 'X'
  
  out = ROOT.TH1D(name, f"#alpha#timesKDE(x);{x_axis_title};Entries", nb, xlo, xhi)
  out.SetDirectory(0)
  for i in range(1, nb + 1):
    x = ref_hist.GetBinCenter(i)
    out.SetBinContent(i, alpha * kde_func.Eval(x))
  return out


def histogram_distribution_stats(hist: ROOT.TH1) -> DistributionStats:
  """Weighted mean, std, and integral using bin centers and contents."""
  integral = 0.0
  mean_num = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    weight = hist.GetBinContent(i)
    if weight <= 0:
      continue
    x = hist.GetBinCenter(i)
    integral += weight
    mean_num += weight * x

  if integral <= 0:
    return DistributionStats(0.0, 0.0, 0.0)

  mean = mean_num / integral
  var_num = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    weight = hist.GetBinContent(i)
    if weight <= 0:
      continue
    x = hist.GetBinCenter(i)
    diff = x - mean
    var_num += weight * diff * diff

  std = math.sqrt(var_num / integral)
  return DistributionStats(mean=mean, std=std, integral=integral)


def _format_stats_row(label: str, stats: DistributionStats, *, width: int = 14) -> str:
  return (
    f"{label:<16}"
    f"{stats.mean:>{width}.6g}"
    f"{stats.std:>{width}.6g}"
    f"{stats.integral:>{width}.6g}"
  )


def print_distribution_stats_table(
  hist_stats: DistributionStats,
  kde_stats: DistributionStats,
  *,
  hist_label: str = "histogram",
  template_label: str = "KDE template",
  unit: str = "cm",
) -> None:
  """Print mean, std, and integral for histogram vs scaled KDE template."""
  delta = DistributionStats(
    mean=kde_stats.mean - hist_stats.mean,
    std=kde_stats.std - hist_stats.std,
    integral=kde_stats.integral - hist_stats.integral,
  )

  col = 14
  print(f"\nDistribution statistics ({unit}):")
  print(
    f"{'':16}{'mean':>{col}}{'std':>{col}}{'integral':>{col}}"
  )
  print(_format_stats_row(hist_label, hist_stats, width=col))
  print(_format_stats_row(template_label, kde_stats, width=col))
  print(_format_stats_row("difference", delta, width=col))


def stats_meta_string(hist_stats: DistributionStats, kde_stats: DistributionStats) -> str:
  return (
    f"hist_mean={hist_stats.mean};hist_std={hist_stats.std};"
    f"hist_integral={hist_stats.integral};"
    f"kde_mean={kde_stats.mean};kde_std={kde_stats.std};"
    f"kde_integral={kde_stats.integral}"
  )


def save_scan_results(
  outfile: str,
  target: ROOT.TH1,
  scan: List[RhoScanPoint],
  g_chi2: ROOT.TGraph,
  g_reduced: ROOT.TGraph,
  g_alpha: ROOT.TGraph,
  best: RhoScanPoint,
  kde_shape: ROOT.TF1,
  template_hist: ROOT.TH1,
  ndf: int,
  *,
  bandwidth: float,
  use_linear_combo: bool = USE_LINEAR_COMBO,
  hist_stats: DistributionStats | None = None,
  kde_stats: DistributionStats | None = None,
) -> None:
  """Write target, scan graphs, best-fit KDE, and metadata."""
  fout = ROOT.TFile.Open(outfile, "RECREATE")
  if not fout or fout.IsZombie():
    raise OSError(f"cannot create output file: {outfile}")

  hist_out = target.Clone("target_hist")
  hist_out.SetDirectory(fout)
  hist_out.Write()

  fout.cd()
  g_chi2.Write()
  g_reduced.Write()
  g_alpha.Write()

  shape = kde_shape.Clone("kde_shape")
  fout.cd()
  shape.Write()

  template_hist.SetDirectory(fout)
  template_hist.Write()

  reduced = best.reduced_chi2
  meta_parts = [
    f"linear_combo={int(use_linear_combo)}",
    f"rho={best.rho}",
    f"bandwidth={bandwidth}",
    f"alpha={best.alpha}",
    f"chi2={best.chi2}",
    f"ndf={ndf}",
    f"reduced_chi2={reduced}",
    f"rho_scan_min={RHO_SCAN_MIN}",
    f"rho_scan_max={RHO_SCAN_MAX}",
    f"rho_scan_step={RHO_SCAN_STEP}",
  ]
  if use_linear_combo:
    meta_parts.insert(2, f"mix={best.mix}")
  meta = ROOT.TNamed("fit_meta", ";".join(meta_parts))
  fout.cd()
  meta.Write()

  if hist_stats is not None and kde_stats is not None:
    stats_meta = ROOT.TNamed("stats_meta", stats_meta_string(hist_stats, kde_stats))
    fout.cd()
    stats_meta.Write()

  fout.Write()
  fout.Close()


def _alpha_y_on_pad(
  alpha: float,
  alpha_lo: float,
  alpha_hi: float,
  pad_lo: float,
  pad_hi: float,
) -> float:
  """Map alpha value to left-pad y coordinate for overlay with chi2/ndf."""
  span = alpha_hi - alpha_lo
  if span <= 0:
    return pad_lo
  frac = (alpha - alpha_lo) / span
  return pad_lo + frac * (pad_hi - pad_lo)


def plot_chi2_vs_bandwidth(
  scan: List[RhoScanPoint],
  g_reduced: ROOT.TGraph,
  ndf: int,
  outfile: str,
) -> None:
  """Save PNG of chi2/ndf vs rho with alpha on a right-hand axis."""
  canvas = ROOT.TCanvas("c_chi2_scan", "chi2 vs bandwidth", 900, 650)
  canvas.SetGrid()
  canvas.SetLeftMargin(0.12)
  canvas.SetRightMargin(0.14)
  canvas.SetBottomMargin(0.12)

  g_reduced.SetLineColor(ROOT.kBlue + 1)
  g_reduced.SetLineWidth(2)
  g_reduced.SetMarkerColor(ROOT.kBlue + 1)
  g_reduced.SetMarkerStyle(20)
  g_reduced.SetMarkerSize(0.7)
  g_reduced.GetYaxis().SetTitleOffset(1.35)
  g_reduced.Draw("ALP")
  canvas.Update()

  pad = canvas.GetPad(0)
  xmax = pad.GetUxmax()
  ymin = pad.GetUymin()
  ymax = pad.GetUymax()

  alphas = [pt.alpha for pt in scan]
  alpha_lo = min(alphas)
  alpha_hi = max(alphas)
  alpha_pad = 0.02 * (alpha_hi - alpha_lo) if alpha_hi > alpha_lo else 0.01 * max(abs(alpha_hi), 1.0)
  alpha_lo -= alpha_pad
  alpha_hi += alpha_pad

  # Vertical axis on the right: same x at both endpoints (not xmin→xmax diagonal).
  axis = ROOT.TGaxis(
    xmax,
    ymin,
    xmax,
    ymax,
    alpha_lo,
    alpha_hi,
    505,
    "-L",
  )
  axis.SetLineColor(ROOT.kRed + 1)
  axis.SetLabelColor(ROOT.kRed + 1)
  axis.SetTitleColor(ROOT.kRed + 1)
  axis.SetTitle("#alpha")
  axis.SetTitleOffset(1.15)
  axis.Draw()

  g_alpha = ROOT.TGraph(len(scan))
  for i, pt in enumerate(scan):
    y_pad = _alpha_y_on_pad(pt.alpha, alpha_lo, alpha_hi, ymin, ymax)
    g_alpha.SetPoint(i, pt.rho, y_pad)
  g_alpha.SetLineColor(ROOT.kRed + 1)
  g_alpha.SetMarkerColor(ROOT.kRed + 1)
  g_alpha.SetMarkerStyle(24)
  g_alpha.SetMarkerSize(0.7)
  g_alpha.SetLineWidth(2)
  g_alpha.Draw("LP SAME")

  leg = ROOT.TLegend(0.52, 0.68, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.03)
  leg.AddEntry(g_reduced, "#chi^{2}/ndf (left)", "lp")
  leg.AddEntry(g_alpha, "#alpha (right)", "lp")
  leg.Draw()

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.03)
  latex.DrawLatex(
    0.14,
    0.92,
    f"#rho #in [{RHO_SCAN_MIN}, {RHO_SCAN_MAX}], step {RHO_SCAN_STEP}; ndf = {ndf}",
  )

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(
    description="Weighted fixed TKDE bandwidth scan with optional mirror blend.",
  )
  parser.add_argument(
    "--linear-combo",
    action=argparse.BooleanOptionalAction,
    default=USE_LINEAR_COMBO,
    help=(
      "optimize a linear blend of unmirrored and mirrored KDEs "
      "(default: %(default)s)"
    ),
  )
  return parser.parse_args()


def main() -> int:
  args = parse_args()
  use_linear_combo = args.linear_combo

  print(f"Reading {TARGET_HIST_NAME!r} from {INPUT_ROOT_FILE}")
  target = open_input_histogram(INPUT_ROOT_FILE, TARGET_HIST_NAME, PROJECTION_AXIS)

  #target = drop_edge_bins(target)
  data, weights = histogram_to_weighted_sample(target)
  x_min = target.GetXaxis().GetXmin()
  x_max = target.GetXaxis().GetXmax()

  nbins = count_chi2_bins(target)
  ndf = degrees_of_freedom(nbins, use_linear_combo=use_linear_combo)

  ctx = KdeFitContext(
    data=data,
    weights=weights,
    x_min=x_min,
    x_max=x_max,
    target=target,
    options=TKDE_OPTIONS,
    ndf=ndf,
  )

  mode = "linear unmirrored/mirrored blend" if use_linear_combo else "single KDE"
  print(f"Weighted TKDE: {len(data)} points on [{x_min}, {x_max}]")
  print(f"Fit mode: {mode}")
  print(f"Options: {TKDE_OPTIONS!r}")
  print(
    f"Bandwidth scan: rho from {RHO_SCAN_MIN} to {RHO_SCAN_MAX} "
    f"step {RHO_SCAN_STEP}  ndf = {ndf}"
  )

  scan = scan_rho_bandwidths(ctx, use_linear_combo=use_linear_combo)

  # print rho for which the reduced chi2 is 1/100 of 1
  for pt in scan:
    if pt.reduced_chi2 < 1.001 and pt.reduced_chi2 > 0.999:
      print("rho: ", pt.rho,"-----", "reduced chi2: ", pt.reduced_chi2)


  best = best_scan_point(scan)
  summary = (
    f"\nMinimum chi2/ndf at rho = {best.rho:.6g}  "
    f"chi2 = {best.chi2:.6g}  chi2/ndf = {best.reduced_chi2:.6g}"
  )
  if use_linear_combo:
    summary = (
      f"\nMinimum chi2/ndf at rho = {best.rho:.6g}  "
      f"mix(unmirrored) = {best.mix:.6g}  "
      f"chi2 = {best.chi2:.6g}  chi2/ndf = {best.reduced_chi2:.6g}"
    )
  print(summary)

  g_chi2, g_reduced, g_alpha = scan_to_graphs(scan)

  kde_unmirrored, kde_mirrored, raw_shape, fit_fn = scaled_kde_tf1(
    best.rho,
    best.alpha,
    ctx,
    mix=best.mix,
    use_linear_combo=use_linear_combo,
  )
  template = kde_template_histogram(raw_shape, target, best.alpha)
  hist_stats = histogram_distribution_stats(target)
  kde_stats = histogram_distribution_stats(template)
  print_distribution_stats_table(hist_stats, kde_stats)

  # print_chi2_contributions(target, raw_shape, best.alpha, exclude_first_last=True)

  if use_linear_combo:
    bandwidth = kde_unmirrored.GetFixedWeight()
    print("bandwidth rho*h_0 (unmirrored): ", bandwidth)
    print("bandwidth rho*h_0 (mirrored): ", kde_mirrored.GetFixedWeight())
  else:
    bandwidth = kde_mirrored.GetFixedWeight()
    print("bandwidth rho*h_0: ", bandwidth)

  out_dir = os.path.dirname(os.path.abspath(OUTPUT_ROOT_FILE))
  os.makedirs(out_dir, exist_ok=True)
  save_scan_results(
    OUTPUT_ROOT_FILE,
    target,
    scan,
    g_chi2,
    g_reduced,
    g_alpha,
    best,
    raw_shape,
    template,
    ndf,
    bandwidth=bandwidth,
    use_linear_combo=use_linear_combo,
    hist_stats=hist_stats,
    kde_stats=kde_stats,
  )
  print(f"Wrote scan and best fit to {OUTPUT_ROOT_FILE}")

  # plot_chi2_vs_bandwidth(scan, g_reduced, ndf, CHI2_SCAN_PNG)

  del kde_unmirrored, kde_mirrored, fit_fn
  return 0


if __name__ == "__main__":
  sys.exit(main())
