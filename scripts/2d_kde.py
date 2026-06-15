#!/usr/bin/env python3
"""
2d_kde.py — Weighted fixed-bandwidth 2D KDE on nominalxyposMM1.

Builds mirrored and unmirrored RooNDKeysPdf models, optimizes their linear blend
and profiled amplitude against the target TH2D, and writes scaled templates to ROOT.
"""

import argparse
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
)
TARGET_HIST_NAME = "nominalxyposMM1"

OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "2d_kde.root"
)

# Bandwidth multiplier for RooNDKeysPdf. Values ~2–3× bin spacing create visible
# lattice ripples on the 8 cm histogram grid; ~1.4 is much smoother.
RHO_SCAN_MIN = 2
RHO_SCAN_MAX = 2
RHO_SCAN_STEP = 0.1

USE_LINEAR_COMBO = False

MIRROR_NO = ROOT.RooNDKeysPdf.NoMirror
MIRROR_BOTH = ROOT.RooNDKeysPdf.MirrorBoth


@dataclass
class RhoScanPoint:
  rho: float
  chi2: float
  reduced_chi2: float
  alpha: float
  mix: float


@dataclass
class KdeFitContext:
  target: ROOT.TH2
  x_var: ROOT.RooRealVar
  y_var: ROOT.RooRealVar
  argset: ROOT.RooArgSet
  dataset: ROOT.RooDataSet
  ndf: int


@dataclass
class Th2Stats:
  integral: float
  mean_x: float
  mean_y: float


def rho_scan_values(
  rho_min: float = RHO_SCAN_MIN,
  rho_max: float = RHO_SCAN_MAX,
  step: float = RHO_SCAN_STEP,
) -> Iterator[float]:
  rho = rho_min
  eps = 1e-12 * max(abs(rho_max), 1.0)
  while rho <= rho_max + eps:
    yield rho
    rho += step


def open_input_histogram(filepath: str, hist_name: str) -> ROOT.TH2:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  obj = tfile.Get(hist_name)
  if not obj:
    tfile.Close()
    raise KeyError(f"object {hist_name!r} not found in {filepath}")

  if not obj.InheritsFrom("TH2"):
    tfile.Close()
    raise TypeError(f"{hist_name!r} is not TH2")

  hist = obj
  hist.SetDirectory(0)
  tfile.Close()
  return hist


def histogram_to_weighted_dataset(hist: ROOT.TH2) -> Tuple[
  ROOT.RooRealVar,
  ROOT.RooRealVar,
  ROOT.RooArgSet,
  ROOT.RooDataSet,
]:
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  ylo = hist.GetYaxis().GetXmin()
  yhi = hist.GetYaxis().GetXmax()

  x_var = ROOT.RooRealVar("x", "x [cm]", xlo, xhi)
  y_var = ROOT.RooRealVar("y", "y [cm]", ylo, yhi)
  argset = ROOT.RooArgSet(x_var, y_var)
  w_var = ROOT.RooRealVar("w", "weight", 0.0, 1.0e20)
  dataset = ROOT.RooDataSet(
    "weighted_points",
    "weighted_points",
    argset,
    ROOT.RooFit.WeightVar(w_var),
  )

  n_added = 0
  for ix in range(1, hist.GetNbinsX() + 1):
    x_var.setVal(hist.GetXaxis().GetBinCenter(ix))
    for iy in range(1, hist.GetNbinsY() + 1):
      content = hist.GetBinContent(ix, iy)
      if content <= 0:
        continue
      y_var.setVal(hist.GetYaxis().GetBinCenter(iy))
      w_var.setVal(content)
      dataset.add(argset, content)
      n_added += 1

  if n_added == 0:
    raise ValueError("histogram has no positive bin content; cannot build weighted KDE")

  return x_var, y_var, argset, dataset


def make_ndkeys_pdf(
  name: str,
  ctx: KdeFitContext,
  *,
  mirror: int,
  rho: float,
) -> ROOT.RooNDKeysPdf:
  return ROOT.RooNDKeysPdf(
    name,
    name,
    ctx.argset,
    ctx.dataset,
    mirror,
    float(rho),
  )


def th2_pdf_shape(
  pdf: ROOT.RooNDKeysPdf,
  ref_hist: ROOT.TH2,
  ctx: KdeFitContext,
  name: str,
) -> ROOT.TH2D:
  """Evaluate an unnormalized RooNDKeysPdf at ref_hist bin centers."""
  nbx = ref_hist.GetNbinsX()
  nby = ref_hist.GetNbinsY()
  out = ROOT.TH2D(
    name,
    name,
    nbx,
    ref_hist.GetXaxis().GetXmin(),
    ref_hist.GetXaxis().GetXmax(),
    nby,
    ref_hist.GetYaxis().GetXmin(),
    ref_hist.GetYaxis().GetXmax(),
  )
  out.SetDirectory(0)

  for ix in range(1, nbx + 1):
    ctx.x_var.setVal(ref_hist.GetXaxis().GetBinCenter(ix))
    for iy in range(1, nby + 1):
      ctx.y_var.setVal(ref_hist.GetYaxis().GetBinCenter(iy))
      out.SetBinContent(ix, iy, pdf.getVal(ctx.argset))

  out._hold_pdf = pdf
  return out


def combined_th2_shape(
  unmirrored: ROOT.TH2,
  mirrored: ROOT.TH2,
  mix: float,
  name: str = "kde_combined_shape",
  title: str = "KDE shape",
) -> ROOT.TH2D:
  out = unmirrored.Clone(name)
  out.SetDirectory(0)
  out.SetTitle(title)
  nbx = out.GetNbinsX()
  nby = out.GetNbinsY()
  for ix in range(1, nbx + 1):
    for iy in range(1, nby + 1):
      u = unmirrored.GetBinContent(ix, iy)
      m = mirrored.GetBinContent(ix, iy)
      out.SetBinContent(ix, iy, mix * u + (1.0 - mix) * m)
  return out


def optimal_alpha(shape: ROOT.TH2, hist: ROOT.TH2) -> float:
  num, den = 0.0, 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      if err <= 0:
        continue
      observed = hist.GetBinContent(ix, iy)
      model_shape = shape.GetBinContent(ix, iy)
      w = 1.0 / (err * err)
      num += w * observed * model_shape
      den += w * model_shape * model_shape
  if den <= 0:
    return 1.0
  return num / den


def chi_squared_vs_hist(
  shape: ROOT.TH2,
  hist: ROOT.TH2,
  alpha: float,
) -> float:
  chi2 = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      if err <= 0:
        continue
      observed = hist.GetBinContent(ix, iy)
      expected = alpha * shape.GetBinContent(ix, iy)
      diff = observed - expected
      chi2 += (diff * diff) / (err * err)
  return chi2


def count_chi2_bins(hist: ROOT.TH2) -> int:
  n = 0
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      if hist.GetBinError(ix, iy) > 0:
        n += 1
  return n


def fit_nparams(use_linear_combo: bool = USE_LINEAR_COMBO) -> int:
  return 2 if use_linear_combo else 1


def degrees_of_freedom(
  nbins: int,
  *,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> int:
  return max(nbins - fit_nparams(use_linear_combo), 1)


def ndkeys_pair_at_rho(
  rho: float,
  ctx: KdeFitContext,
) -> Tuple[ROOT.RooNDKeysPdf, ROOT.TH2D, ROOT.RooNDKeysPdf, ROOT.TH2D]:
  pdf_unmirrored = make_ndkeys_pdf(
    f"kde_unmirrored_rho{rho:.4g}",
    ctx,
    mirror=MIRROR_NO,
    rho=rho,
  )
  pdf_mirrored = make_ndkeys_pdf(
    f"kde_mirrored_rho{rho:.4g}",
    ctx,
    mirror=MIRROR_BOTH,
    rho=rho,
  )
  th2_unmirrored = th2_pdf_shape(
    pdf_unmirrored,
    ctx.target,
    ctx,
    f"kde_unmirrored_shape_rho{rho:.4g}",
  )
  th2_mirrored = th2_pdf_shape(
    pdf_mirrored,
    ctx.target,
    ctx,
    f"kde_mirrored_shape_rho{rho:.4g}",
  )
  return pdf_unmirrored, th2_unmirrored, pdf_mirrored, th2_mirrored


def evaluate_at_mix(
  mix: float,
  th2_unmirrored: ROOT.TH2,
  th2_mirrored: ROOT.TH2,
  target: ROOT.TH2,
  ndf: int,
) -> Tuple[float, float, float]:
  shape = combined_th2_shape(th2_unmirrored, th2_mirrored, mix)
  alpha = optimal_alpha(shape, target)
  chi2 = chi_squared_vs_hist(shape, target, alpha)
  reduced = chi2 / max(ndf, 1)
  return chi2, reduced, alpha


def _golden_section_minimize(
  func,
  lo: float,
  hi: float,
  *,
  tol: float = 1e-4,
  max_iter: int = 80,
) -> float:
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
  th2_unmirrored: ROOT.TH2,
  th2_mirrored: ROOT.TH2,
  target: ROOT.TH2,
  ndf: int,
) -> Tuple[float, float, float, float]:
  cache: dict[float, Tuple[float, float, float]] = {}

  def loss(mix: float) -> float:
    if mix not in cache:
      cache[mix] = evaluate_at_mix(mix, th2_unmirrored, th2_mirrored, target, ndf)
    _chi2, reduced, _alpha = cache[mix]
    if not math.isfinite(reduced):
      return math.inf
    return reduced

  mix_opt = _golden_section_minimize(loss, 0.0, 1.0)
  chi2, reduced, alpha = evaluate_at_mix(
    mix_opt,
    th2_unmirrored,
    th2_mirrored,
    target,
    ndf,
  )
  return mix_opt, alpha, chi2, reduced


def evaluate_at_rho(
  rho: float,
  ctx: KdeFitContext,
  *,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> Tuple[float, float, float, float, ROOT.TH2D | None]:
  if not use_linear_combo:
    pdf = make_ndkeys_pdf("kde_single", ctx, mirror=MIRROR_BOTH, rho=rho)
    shape = th2_pdf_shape(pdf, ctx.target, ctx, f"kde_shape_rho{rho:.4g}")
    alpha = optimal_alpha(shape, ctx.target)
    chi2 = chi_squared_vs_hist(shape, ctx.target, alpha)
    reduced = chi2 / max(ctx.ndf, 1)
    return chi2, reduced, alpha, 0.0, shape

  _pdf_u, th2_u, _pdf_m, th2_m = ndkeys_pair_at_rho(rho, ctx)
  mix, alpha, chi2, reduced = optimize_mix_weight(
    th2_u,
    th2_m,
    ctx.target,
    ctx.ndf,
  )
  shape = combined_th2_shape(th2_u, th2_m, mix, name=f"kde_shape_rho{rho:.4g}")
  shape._hold_unmirrored = th2_u
  shape._hold_mirrored = th2_m
  return chi2, reduced, alpha, mix, shape


def scan_rho_bandwidths(
  ctx: KdeFitContext,
  *,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> List[RhoScanPoint]:
  results: List[RhoScanPoint] = []
  for rho in rho_scan_values():
    chi2, reduced, alpha, mix, _shape = evaluate_at_rho(
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
  return min(scan, key=lambda p: p.reduced_chi2)


def scaled_kde_th2(
  rho: float,
  alpha: float,
  ctx: KdeFitContext,
  *,
  mix: float = 1.0,
  use_linear_combo: bool = USE_LINEAR_COMBO,
) -> Tuple[ROOT.TH2D, ROOT.TH2D]:
  """Return (shape, alpha * shape) TH2 templates at the best scan point."""
  if not use_linear_combo:
    pdf = make_ndkeys_pdf("kde_single_best", ctx, mirror=MIRROR_BOTH, rho=rho)
    shape = th2_pdf_shape(pdf, ctx.target, ctx, "kde_shape")
    template = shape.Clone("kde_template")
    template.SetDirectory(0)
    template.Scale(alpha)
    template.SetTitle("#alpha#timesKDE(x,y);x [cm];y [cm]")
    return shape, template

  pdf_u, th2_u, pdf_m, th2_m = ndkeys_pair_at_rho(rho, ctx)
  shape = combined_th2_shape(
    th2_u,
    th2_m,
    mix,
    name="kde_shape",
    title="mirrored/unmirrored KDE shape",
  )
  shape._hold_pdf_unmirrored = pdf_u
  shape._hold_pdf_mirrored = pdf_m
  shape._hold_unmirrored = th2_u
  shape._hold_mirrored = th2_m

  template = shape.Clone("kde_template")
  template.SetDirectory(0)
  template.Scale(alpha)
  template.SetTitle("#alpha#timesKDE(x,y);x [cm];y [cm]")
  return shape, template


def th2_distribution_stats(hist: ROOT.TH2) -> Th2Stats:
  integral = 0.0
  mean_x_num = 0.0
  mean_y_num = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    x = hist.GetXaxis().GetBinCenter(ix)
    for iy in range(1, hist.GetNbinsY() + 1):
      weight = hist.GetBinContent(ix, iy)
      if weight <= 0:
        continue
      y = hist.GetYaxis().GetBinCenter(iy)
      integral += weight
      mean_x_num += weight * x
      mean_y_num += weight * y

  if integral <= 0:
    return Th2Stats(0.0, 0.0, 0.0)

  return Th2Stats(
    integral=integral,
    mean_x=mean_x_num / integral,
    mean_y=mean_y_num / integral,
  )


def print_distribution_stats_table(
  hist_stats: Th2Stats,
  kde_stats: Th2Stats,
) -> None:
  print("\nDistribution statistics:")
  print(f"{'':16}{'integral':>16}{'mean_x':>16}{'mean_y':>16}")
  print(
    f"{'histogram':<16}"
    f"{hist_stats.integral:16.6g}"
    f"{hist_stats.mean_x:16.6g}"
    f"{hist_stats.mean_y:16.6g}"
  )
  print(
    f"{'KDE template':<16}"
    f"{kde_stats.integral:16.6g}"
    f"{kde_stats.mean_x:16.6g}"
    f"{kde_stats.mean_y:16.6g}"
  )
  print(
    f"{'difference':<16}"
    f"{kde_stats.integral - hist_stats.integral:16.6g}"
    f"{kde_stats.mean_x - hist_stats.mean_x:16.6g}"
    f"{kde_stats.mean_y - hist_stats.mean_y:16.6g}"
  )


def stats_meta_string(hist_stats: Th2Stats, kde_stats: Th2Stats) -> str:
  return (
    f"hist_integral={hist_stats.integral};"
    f"hist_mean_x={hist_stats.mean_x};hist_mean_y={hist_stats.mean_y};"
    f"kde_integral={kde_stats.integral};"
    f"kde_mean_x={kde_stats.mean_x};kde_mean_y={kde_stats.mean_y}"
  )


def save_results(
  outfile: str,
  target: ROOT.TH2,
  scan: List[RhoScanPoint],
  g_chi2: ROOT.TGraph,
  g_reduced: ROOT.TGraph,
  g_alpha: ROOT.TGraph,
  best: RhoScanPoint,
  kde_shape: ROOT.TH2,
  template_hist: ROOT.TH2,
  ndf: int,
  *,
  use_linear_combo: bool = USE_LINEAR_COMBO,
  hist_stats: Th2Stats | None = None,
  kde_stats: Th2Stats | None = None,
) -> None:
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
  shape.SetDirectory(fout)
  shape.Write()

  template_hist.SetDirectory(fout)
  template_hist.Write()

  meta_parts = [
    f"linear_combo={int(use_linear_combo)}",
    f"rho={best.rho}",
    f"alpha={best.alpha}",
    f"chi2={best.chi2}",
    f"ndf={ndf}",
    f"reduced_chi2={best.reduced_chi2}",
    f"rho_scan_min={RHO_SCAN_MIN}",
    f"rho_scan_max={RHO_SCAN_MAX}",
    f"rho_scan_step={RHO_SCAN_STEP}",
    "pdf=RooNDKeysPdf",
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


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(
    description="Weighted fixed 2D RooNDKeysPdf with optional mirror blend.",
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
  target = open_input_histogram(INPUT_ROOT_FILE, TARGET_HIST_NAME)
  x_var, y_var, argset, dataset = histogram_to_weighted_dataset(target)

  nbins = count_chi2_bins(target)
  ndf = degrees_of_freedom(nbins, use_linear_combo=use_linear_combo)

  ctx = KdeFitContext(
    target=target,
    x_var=x_var,
    y_var=y_var,
    argset=argset,
    dataset=dataset,
    ndf=ndf,
  )

  mode = "linear unmirrored/mirrored blend" if use_linear_combo else "single KDE"
  print(
    f"Weighted RooNDKeysPdf: {dataset.numEntries()} points on "
    f"x=[{target.GetXaxis().GetXmin()}, {target.GetXaxis().GetXmax()}], "
    f"y=[{target.GetYaxis().GetXmin()}, {target.GetYaxis().GetXmax()}]"
  )
  print(f"Fit mode: {mode}")
  print(
    f"Bandwidth scan: rho from {RHO_SCAN_MIN} to {RHO_SCAN_MAX} "
    f"step {RHO_SCAN_STEP}  ndf = {ndf}"
  )

  scan = scan_rho_bandwidths(ctx, use_linear_combo=use_linear_combo)
  best = best_scan_point(scan)

  if use_linear_combo:
    print(
      f"\nMinimum chi2/ndf at rho = {best.rho:.6g}  "
      f"mix(unmirrored) = {best.mix:.6g}  "
      f"alpha = {best.alpha:.6g}  "
      f"chi2 = {best.chi2:.6g}  chi2/ndf = {best.reduced_chi2:.6g}"
    )
  else:
    print(
      f"\nMinimum chi2/ndf at rho = {best.rho:.6g}  "
      f"alpha = {best.alpha:.6g}  "
      f"chi2 = {best.chi2:.6g}  chi2/ndf = {best.reduced_chi2:.6g}"
    )

  g_chi2, g_reduced, g_alpha = scan_to_graphs(scan)
  kde_shape, template = scaled_kde_th2(
    best.rho,
    best.alpha,
    ctx,
    mix=best.mix,
    use_linear_combo=use_linear_combo,
  )

  hist_stats = th2_distribution_stats(target)
  kde_stats = th2_distribution_stats(template)
  print_distribution_stats_table(hist_stats, kde_stats)

  out_dir = os.path.dirname(os.path.abspath(OUTPUT_ROOT_FILE))
  os.makedirs(out_dir, exist_ok=True)
  save_results(
    OUTPUT_ROOT_FILE,
    target,
    scan,
    g_chi2,
    g_reduced,
    g_alpha,
    best,
    kde_shape,
    template,
    ndf,
    use_linear_combo=use_linear_combo,
    hist_stats=hist_stats,
    kde_stats=kde_stats,
  )
  print(f"Wrote best fit to {OUTPUT_ROOT_FILE}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
