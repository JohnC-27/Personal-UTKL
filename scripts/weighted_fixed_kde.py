#!/usr/bin/env python3
"""
weighted_adaptive_kde.py — Weighted adaptive ROOT.TKDE bandwidth scan.

Scans bandwidth scale rho from RHO_SCAN_MIN to RHO_SCAN_MAX, profiles alpha at each
point, records chi-squared vs the target TH1, writes results to ROOT, and plots
chi2/ndf vs bandwidth.
"""
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
)
TARGET_HIST_NAME = "nominalxyposMM1"  # TH2D; projected to X below
PROJECTION_AXIS = "x"  # "x" or "y" for TH2 projections; ignored for TH1

OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "weighted_fixed_kde.root"
)
CHI2_SCAN_PNG = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "weighted_fixed_kde_chi2_scan.png"
)

TKDE_OPTIONS = "KernelType:Gaussian;Iteration:Fixed;Mirror:noMirror"

NPAR = 1  # rho; alpha profiled (linear WLS)

RHO_SCAN_MIN = 0.0983
RHO_SCAN_MAX = 0.1
RHO_SCAN_STEP = 0.00001

# optimal rho for chi2/ndf = 1
# mm1x: 0.0983   y: 0.10486 or 0.10487
# mm2x: 




@dataclass
class RhoScanPoint:
  rho: float
  chi2: float
  reduced_chi2: float
  alpha: float


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


def chi2_contributions(hist: ROOT.TH1, kde_func: ROOT.TF1, alpha: float) -> list:
  """Makes a list of the contributions to the chi2 error"""
  contributions = []
  for i in range(1, count_chi2_bins(hist)):
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    observed = hist.GetBinContent(i)
    expected = alpha * kde_func.Eval(hist.GetBinCenter(i))
    diff = observed - expected
    contributions[i] = (diff * diff) / (err * err)
  return contributions


def degrees_of_freedom(nbins: int, npar: int = NPAR) -> int:
  """ndf = nbins - nparams."""
  return max(nbins - npar, 1)


def evaluate_at_rho(rho: float, ctx: KdeFitContext) -> Tuple[float, float, float]:
  """Return (chi2, reduced chi2, profiled alpha) at fixed rho."""
  kde = make_tkde(rho, ctx)
  kde_func = kde.GetFunction()
  if kde_func is None:
    return math.inf, math.inf, 1.0
  kde_func.SetNpx(ctx.npx)
  alpha = optimal_alpha(kde_func, ctx.target)
  chi2 = chi_squared_vs_hist(kde_func, ctx.target, alpha)
  reduced = chi2 / max(ctx.ndf, 1)
  return chi2, reduced, alpha


def scan_rho_bandwidths(ctx: KdeFitContext) -> List[RhoScanPoint]:
  """Evaluate chi2 at each rho on the scan grid."""
  results: List[RhoScanPoint] = []
  for rho in rho_scan_values():
    chi2, reduced, alpha = evaluate_at_rho(rho, ctx)
    if not math.isfinite(chi2):
      print(f"  rho={rho:.4g}: KDE failed, skipped")
      continue
    results.append(RhoScanPoint(rho=rho, chi2=chi2, reduced_chi2=reduced, alpha=alpha))
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
) -> Tuple[ROOT.TKDE, ROOT.TF1]:
  """Build weighted TKDE and TF1 wrapper f(x) = alpha * kde(x)."""
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
  return kde, fit_fn


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
  meta = ROOT.TNamed(
    "fit_meta",
    (
      f"rho={best.rho};alpha={best.alpha};chi2={best.chi2};ndf={ndf};"
      f"reduced_chi2={reduced};rho_scan_min={RHO_SCAN_MIN};"
      f"rho_scan_max={RHO_SCAN_MAX};rho_scan_step={RHO_SCAN_STEP}"
    ),
  )
  fout.cd()
  meta.Write()

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
  print(
    f"Bandwidth scan: rho from {RHO_SCAN_MIN} to {RHO_SCAN_MAX} "
    f"step {RHO_SCAN_STEP}  ndf = {ndf}"
  )

  scan = scan_rho_bandwidths(ctx)

  # print rho for which the reduced chi2 is 1/100 of 1
  for pt in scan:
    if pt.reduced_chi2 < 1.001 and pt.reduced_chi2 > 0.999:
      print("rho: ", pt.rho,"-----", "reduced chi2: ", pt.reduced_chi2)


  best = best_scan_point(scan)
  print(
    f"\nMinimum chi2/ndf at rho = {best.rho:.6g}  "
    f"chi2 = {best.chi2:.6g}  chi2/ndf = {best.reduced_chi2:.6g}"
  )

  g_chi2, g_reduced, g_alpha = scan_to_graphs(scan)

  kde, fit_fn = scaled_kde_tf1(best.rho, best.alpha, ctx)
  raw_shape = fit_fn._hold_raw
  template = kde_template_histogram(raw_shape, target, best.alpha)

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
  )
  print(f"Wrote scan and best fit to {OUTPUT_ROOT_FILE}")

  plot_chi2_vs_bandwidth(scan, g_reduced, ndf, CHI2_SCAN_PNG)
  del kde, fit_fn
  return 0


if __name__ == "__main__":
  sys.exit(main())
