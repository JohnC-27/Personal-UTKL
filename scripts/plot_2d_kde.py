#!/usr/bin/env python3
"""Plot 2D KDE fit results from 2d_kde.root."""

import array
import json
import os
import sys
import time
from dataclasses import dataclass

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

FIT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "2d_kde.root"
)
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
OUTPUT_OVERLAY = os.path.join(OUTPUT_DIR, "2d_kde_overlay.png")
OUTPUT_PROJECTIONS = os.path.join(OUTPUT_DIR, "2d_kde_projections.png")
OUTPUT_RATIO = os.path.join(OUTPUT_DIR, "2d_kde_ratio.png")

RATIO_Z_PAD = 1.05
RATIO_Z_MIN_HALF_WIDTH = 0.05

# Direct RooNDKeysPdf evaluation grid for surfaces and projections.
KDE_PLOT_BINS = 200

# 2D RooNDKeysPdf(RooArgSet, ...) takes an options string, not the legacy Mirror enum.
NDKEYS_NO_MIRROR = "a"
NDKEYS_MIRROR_BOTH = "am"

DEBUG_LOG_PATH = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), ".cursor", "debug-8940bd.log"
)


@dataclass
class KdeEvalContext:
  x_var: ROOT.RooRealVar
  y_var: ROOT.RooRealVar
  argset: ROOT.RooArgSet
  dataset: ROOT.RooDataSet


@dataclass
class KdeModel:
  ctx: KdeEvalContext
  alpha: float
  mix: float
  use_linear_combo: bool
  pdf_single: ROOT.RooNDKeysPdf | None = None
  pdf_unmirrored: ROOT.RooNDKeysPdf | None = None
  pdf_mirrored: ROOT.RooNDKeysPdf | None = None

  def shape_at(self, x: float, y: float) -> float:
    self.ctx.x_var.setVal(x)
    self.ctx.y_var.setVal(y)
    if self.use_linear_combo:
      u = self.pdf_unmirrored.getVal(self.ctx.argset)
      m = self.pdf_mirrored.getVal(self.ctx.argset)
      return self.mix * u + (1.0 - self.mix) * m
    return self.pdf_single.getVal(self.ctx.argset)

  def scaled_at(self, x: float, y: float) -> float:
    return self.alpha * self.shape_at(x, y)


@dataclass
class Th2Stats:
  integral: float
  mean_x: float
  mean_y: float


def parse_meta(meta: ROOT.TNamed) -> dict[str, float | str | int]:
  out: dict[str, float | str | int] = {}
  for part in meta.GetTitle().split(";"):
    key, val = part.split("=", 1)
    if key == "pdf":
      out[key] = val
      continue
    try:
      num = float(val)
      out[key] = int(num) if num.is_integer() and "." not in val else num
    except ValueError:
      out[key] = val
  return out


def parse_stats_meta(meta: ROOT.TNamed) -> tuple[Th2Stats, Th2Stats]:
  raw = parse_meta(meta)
  hist_stats = Th2Stats(
    integral=float(raw["hist_integral"]),
    mean_x=float(raw["hist_mean_x"]),
    mean_y=float(raw["hist_mean_y"]),
  )
  kde_stats = Th2Stats(
    integral=float(raw["kde_integral"]),
    mean_x=float(raw["kde_mean_x"]),
    mean_y=float(raw["kde_mean_y"]),
  )
  return hist_stats, kde_stats


def histogram_to_weighted_dataset(hist: ROOT.TH2) -> KdeEvalContext:
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  ylo = hist.GetYaxis().GetXmin()
  yhi = hist.GetYaxis().GetXmax()

  x_var = ROOT.RooRealVar("plot_x", "x [cm]", xlo, xhi)
  y_var = ROOT.RooRealVar("plot_y", "y [cm]", ylo, yhi)
  argset = ROOT.RooArgSet(x_var, y_var)
  w_var = ROOT.RooRealVar("plot_w", "weight", 0.0, 1.0e20)
  dataset = ROOT.RooDataSet(
    "plot_weighted_points",
    "plot_weighted_points",
    argset,
    ROOT.RooFit.WeightVar(w_var),
  )

  for ix in range(1, hist.GetNbinsX() + 1):
    x_var.setVal(hist.GetXaxis().GetBinCenter(ix))
    for iy in range(1, hist.GetNbinsY() + 1):
      content = hist.GetBinContent(ix, iy)
      if content <= 0:
        continue
      y_var.setVal(hist.GetYaxis().GetBinCenter(iy))
      w_var.setVal(content)
      dataset.add(argset, content)

  return KdeEvalContext(x_var=x_var, y_var=y_var, argset=argset, dataset=dataset)


def _weighted_mean_std_from_th2(hist: ROOT.TH2) -> tuple[float, float, float, float, float]:
  total_w = 0.0
  sum_x = 0.0
  sum_y = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    x = hist.GetXaxis().GetBinCenter(ix)
    for iy in range(1, hist.GetNbinsY() + 1):
      w = hist.GetBinContent(ix, iy)
      if w <= 0:
        continue
      y = hist.GetYaxis().GetBinCenter(iy)
      total_w += w
      sum_x += w * x
      sum_y += w * y

  if total_w <= 0:
    return 0.0, 0.0, 0.0, 0.0, 0.0

  mean_x = sum_x / total_w
  mean_y = sum_y / total_w

  var_x = 0.0
  var_y = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    x = hist.GetXaxis().GetBinCenter(ix)
    dx2 = (x - mean_x) * (x - mean_x)
    for iy in range(1, hist.GetNbinsY() + 1):
      w = hist.GetBinContent(ix, iy)
      if w <= 0:
        continue
      y = hist.GetYaxis().GetBinCenter(iy)
      dy2 = (y - mean_y) * (y - mean_y)
      var_x += w * dx2
      var_y += w * dy2

  std_x = (var_x / total_w) ** 0.5
  std_y = (var_y / total_w) ** 0.5
  return total_w, mean_x, mean_y, std_x, std_y


def ndkeys_bandwidths(target: ROOT.TH2, rho: float) -> tuple[float, float]:
  """Return effective RooNDKeys bandwidths (rho * h0) in x and y."""
  n_eff, _mx, _my, sigma_x, sigma_y = _weighted_mean_std_from_th2(target)
  if n_eff <= 0:
    return 0.0, 0.0

  d = 2.0
  silverman = (4.0 / (d + 2.0)) ** (1.0 / (d + 4.0))
  n_factor = n_eff ** (-1.0 / (d + 4.0))
  h0_x = silverman * sigma_x * n_factor
  h0_y = silverman * sigma_y * n_factor
  return rho * h0_x, rho * h0_y


def print_ndkeys_bandwidths(target: ROOT.TH2, rho: float) -> None:
  """Print RooNDKeys Silverman-rule h0 and effective rho*h0 widths."""
  n_eff, _mx, _my, sigma_x, sigma_y = _weighted_mean_std_from_th2(target)
  if n_eff <= 0:
    print("RooNDKeys bandwidths: unable to compute (empty histogram).")
    return

  d = 2.0
  silverman = (4.0 / (d + 2.0)) ** (1.0 / (d + 4.0))
  n_factor = n_eff ** (-1.0 / (d + 4.0))
  h0_x = silverman * sigma_x * n_factor
  h0_y = silverman * sigma_y * n_factor
  bw_x, bw_y = ndkeys_bandwidths(target, rho)
  print(
    "RooNDKeys bandwidths (Silverman): "
    f"n_eff={n_eff:.6g}, rho={rho:.6g}, "
    f"h0_x={h0_x:.6g}, h0_y={h0_y:.6g}, "
    f"h_x={bw_x:.6g}, h_y={bw_y:.6g}"
  )


def make_ndkeys_pdf(
  name: str,
  ctx: KdeEvalContext,
  *,
  mirror_options: str,
  rho: float,
) -> ROOT.RooNDKeysPdf:
  return ROOT.RooNDKeysPdf(
    name,
    name,
    ctx.argset,
    ctx.dataset,
    mirror_options,
    float(rho),
  )


def _ndkeys_options(meta: dict[str, float | str | int], key: str, default: str) -> str:
  value = meta.get(key, default)
  return str(value)


def build_kde_model(target: ROOT.TH2, meta: dict[str, float | str | int]) -> KdeModel:
  ctx = histogram_to_weighted_dataset(target)
  rho = float(meta["rho"])
  alpha = float(meta["alpha"])
  use_linear_combo = bool(meta.get("linear_combo", 0))
  mix = float(meta.get("mix", 1.0))
  opt_no_mirror = _ndkeys_options(meta, "ndkeys_no_mirror", NDKEYS_NO_MIRROR)
  opt_mirror = _ndkeys_options(meta, "ndkeys_mirror", NDKEYS_MIRROR_BOTH)

  if use_linear_combo:
    pdf_unmirrored = make_ndkeys_pdf(
      "plot_kde_unmirrored",
      ctx,
      mirror_options=opt_no_mirror,
      rho=rho,
    )
    pdf_mirrored = make_ndkeys_pdf(
      "plot_kde_mirrored",
      ctx,
      mirror_options=opt_mirror,
      rho=rho,
    )
    return KdeModel(
      ctx=ctx,
      alpha=alpha,
      mix=mix,
      use_linear_combo=True,
      pdf_unmirrored=pdf_unmirrored,
      pdf_mirrored=pdf_mirrored,
    )

  pdf_single = make_ndkeys_pdf(
    "plot_kde_single",
    ctx,
    mirror_options=opt_mirror,
    rho=rho,
  )
  return KdeModel(
    ctx=ctx,
    alpha=alpha,
    mix=1.0,
    use_linear_combo=False,
    pdf_single=pdf_single,
  )


def evaluate_kde_th2(
  model: KdeModel,
  *,
  n_bins_x: int,
  n_bins_y: int,
  xlo: float,
  xhi: float,
  ylo: float,
  yhi: float,
  name: str,
  title: str,
  scaled: bool = True,
) -> ROOT.TH2D:
  out = ROOT.TH2D(name, title, n_bins_x, xlo, xhi, n_bins_y, ylo, yhi)
  out.SetDirectory(0)
  out.SetStats(0)
  out._hold_model = model

  value_at = model.scaled_at if scaled else model.shape_at
  for ix in range(1, n_bins_x + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    for iy in range(1, n_bins_y + 1):
      y = out.GetYaxis().GetBinCenter(iy)
      out.SetBinContent(ix, iy, value_at(x, y))

  return out


def evaluate_kde_on_hist_grid(
  model: KdeModel,
  ref_hist: ROOT.TH2,
  name: str,
  *,
  scaled: bool = True,
) -> ROOT.TH2D:
  return evaluate_kde_th2(
    model,
    n_bins_x=ref_hist.GetNbinsX(),
    n_bins_y=ref_hist.GetNbinsY(),
    xlo=ref_hist.GetXaxis().GetXmin(),
    xhi=ref_hist.GetXaxis().GetXmax(),
    ylo=ref_hist.GetYaxis().GetXmin(),
    yhi=ref_hist.GetYaxis().GetXmax(),
    name=name,
    title=ref_hist.GetTitle(),
    scaled=scaled,
  )


def kde_plot_hist(model: KdeModel, name: str) -> ROOT.TH2D:
  target = model.ctx.x_var
  xlo = target.getMin()
  xhi = target.getMax()
  ylo = model.ctx.y_var.getMin()
  yhi = model.ctx.y_var.getMax()
  return evaluate_kde_th2(
    model,
    n_bins_x=KDE_PLOT_BINS,
    n_bins_y=KDE_PLOT_BINS,
    xlo=xlo,
    xhi=xhi,
    ylo=ylo,
    yhi=yhi,
    name=name,
    title="#alpha#timesKDE(x,y);x [cm];y [cm]",
    scaled=True,
  )


def _hist_axis_centers(axis: ROOT.TAxis) -> list[float]:
  return [axis.GetBinCenter(i) for i in range(1, axis.GetNbins() + 1)]


def _axis_centers(lo: float, hi: float, n_bins: int) -> list[float]:
  width = (hi - lo) / n_bins
  return [lo + (i - 0.5) * width for i in range(1, n_bins + 1)]


def kde_projection_x(
  model: KdeModel,
  ref_hist: ROOT.TH2,
  name: str,
  *,
  n_bins: int = KDE_PLOT_BINS,
  n_integrate: int = KDE_PLOT_BINS,
) -> ROOT.TH1D:
  """X marginal from direct PDF evaluation at many x, summed over y."""
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()
  ylo = ref_hist.GetYaxis().GetXmin()
  yhi = ref_hist.GetYaxis().GetXmax()
  out = ROOT.TH1D(name, name, n_bins, xlo, xhi)
  out.SetDirectory(0)
  out.SetStats(0)
  out._hold_model = model

  y_values = _axis_centers(ylo, yhi, n_integrate)
  y_scale = ref_hist.GetNbinsY() / n_integrate
  for ix in range(1, n_bins + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    marginal = sum(model.scaled_at(x, y) for y in y_values) * y_scale
    out.SetBinContent(ix, marginal)

  return out


def kde_projection_y(
  model: KdeModel,
  ref_hist: ROOT.TH2,
  name: str,
  *,
  n_bins: int = KDE_PLOT_BINS,
  n_integrate: int = KDE_PLOT_BINS,
) -> ROOT.TH1D:
  """Y marginal from direct PDF evaluation at many y, summed over x."""
  ylo = ref_hist.GetYaxis().GetXmin()
  yhi = ref_hist.GetYaxis().GetXmax()
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()
  out = ROOT.TH1D(name, name, n_bins, ylo, yhi)
  out.SetDirectory(0)
  out.SetStats(0)
  out._hold_model = model

  x_values = _axis_centers(xlo, xhi, n_integrate)
  x_scale = ref_hist.GetNbinsX() / n_integrate
  for iy in range(1, n_bins + 1):
    y = out.GetXaxis().GetBinCenter(iy)
    marginal = sum(model.scaled_at(x, y) for x in x_values) * x_scale
    out.SetBinContent(iy, marginal)

  return out


def load_fit_objects(filepath: str):
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  target = tfile.Get("target_hist")
  template = tfile.Get("kde_template")
  fit_meta = tfile.Get("fit_meta")
  stats_meta = tfile.Get("stats_meta")
  for name, obj in [
    ("target_hist", target),
    ("kde_template", template),
    ("fit_meta", fit_meta),
    ("stats_meta", stats_meta),
  ]:
    if not obj:
      tfile.Close()
      raise KeyError(f"missing {name!r} in {filepath}")

  target = target.Clone("target_hist_plot")
  template = template.Clone("kde_template_plot")
  target.SetDirectory(0)
  template.SetDirectory(0)
  target.SetStats(0)
  template.SetStats(0)

  meta = parse_meta(fit_meta)
  hist_stats, kde_stats = parse_stats_meta(stats_meta)
  tfile.Close()
  return target, template, meta, hist_stats, kde_stats


def _axis_titles(hist: ROOT.TH2) -> tuple[str, str, str]:
  return (
    hist.GetXaxis().GetTitle() or "x [cm]",
    hist.GetYaxis().GetTitle() or "y [cm]",
    hist.GetZaxis().GetTitle() or "Entries",
  )


def _style_surf3d_axes(hist: ROOT.TH2) -> tuple[str, str, str]:
  xtitle, ytitle, ztitle = _axis_titles(hist)
  for axis in (hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()):
    axis.SetTitle("")
    axis.SetLabelSize(0.022)
    axis.SetNdivisions(505)
  return xtitle, ytitle, ztitle


def _draw_surf3d_axis_titles(xtitle: str, ytitle: str, ztitle: str) -> None:
  latex = ROOT.TLatex()
  latex.SetNDC(True)
  latex.SetTextFont(42)
  latex.SetTextSize(0.032)
  latex.SetTextAlign(22)
  latex.DrawLatex(0.74, 0.05, xtitle)
  latex.DrawLatex(0.26, 0.05, ytitle)
  latex.SetTextAngle(90)
  latex.DrawLatex(0.97, 0.54, ztitle)


def _configure_surf_canvas(pad: ROOT.TPad) -> None:
  pad.SetGrid()
  pad.SetTheta(28)
  pad.SetPhi(60)
  pad.SetLeftMargin(0.07)
  pad.SetRightMargin(0.07)
  pad.SetBottomMargin(0.05)
  pad.SetTopMargin(0.10)


def _style_surf_hist(
  hist: ROOT.TH2,
  *,
  line_color: int,
  line_width: int = 1,
) -> tuple[str, str, str]:
  hist.SetLineColor(line_color)
  hist.SetLineWidth(line_width)
  hist.SetFillStyle(0)
  return _style_surf3d_axes(hist)


def _set_diverging_ratio_palette() -> None:
  """Blue (ratio < 1) -> white (ratio = 1) -> red (ratio > 1)."""
  stops = array.array("d", [0.0, 0.5, 1.0])
  red = array.array("d", [0.0, 1.0, 1.0])
  green = array.array("d", [0.0, 1.0, 0.0])
  blue = array.array("d", [1.0, 1.0, 0.0])
  ROOT.TColor.CreateGradientColorTable(3, stops, red, green, blue, 255)
  ROOT.gStyle.SetNumberContours(255)


def kde_over_data_ratio(kde: ROOT.TH2, data: ROOT.TH2, name: str) -> ROOT.TH2D:
  ratio = kde.Clone(name)
  ratio.SetDirectory(0)
  ratio.Divide(data)
  ratio.SetStats(0)
  return ratio


def _ratio_color_range(
  ratio: ROOT.TH2,
  data: ROOT.TH2,
  *,
  center: float = 1.0,
  pad: float = RATIO_Z_PAD,
  min_half_width: float = RATIO_Z_MIN_HALF_WIDTH,
) -> tuple[float, float]:
  max_dev = 0.0
  for ix in range(1, data.GetNbinsX() + 1):
    for iy in range(1, data.GetNbinsY() + 1):
      if data.GetBinContent(ix, iy) <= 0:
        continue
      max_dev = max(max_dev, abs(ratio.GetBinContent(ix, iy) - center))

  half = max(max_dev * pad, min_half_width)
  return center - half, center + half


def plot_ratio(
  target: ROOT.TH2,
  template: ROOT.TH2,
  outfile: str,
) -> None:
  ratio = kde_over_data_ratio(template, target, "kde_over_data_ratio")
  zmin, zmax = _ratio_color_range(ratio, target)

  xtitle = target.GetXaxis().GetTitle() or "x [cm]"
  ytitle = target.GetYaxis().GetTitle() or "y [cm]"
  ratio.SetTitle("KDE / Data;{};{};KDE / Data".format(xtitle, ytitle))
  ratio.GetZaxis().SetRangeUser(zmin, zmax)
  ratio.GetZaxis().SetTitle("KDE / Data")
  ratio.GetXaxis().SetTitle(xtitle)
  ratio.GetYaxis().SetTitle(ytitle)

  _set_diverging_ratio_palette()

  canvas = ROOT.TCanvas("c_ratio", "KDE / Data ratio", 900, 820)
  canvas.SetRightMargin(0.14)
  canvas.SetLeftMargin(0.12)
  canvas.SetBottomMargin(0.12)
  ratio.Draw("COLZ")

  latex = ROOT.TLatex()
  latex.SetNDC(True)
  latex.SetTextFont(42)
  latex.SetTextSize(0.03)
  latex.SetTextAlign(22)
  latex.DrawLatex(
    0.50,
    0.92,
    "blue: KDE < data   white: KDE = data   red: KDE > data",
  )

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def _draw_stats_and_params(
  pad: ROOT.TPad,
  target: ROOT.TH2,
  hist_stats: Th2Stats,
  kde_stats: Th2Stats,
  meta: dict[str, float | str | int],
) -> None:
  pad.cd()
  pad.SetFillStyle(0)

  delta = Th2Stats(
    integral=kde_stats.integral - hist_stats.integral,
    mean_x=kde_stats.mean_x - hist_stats.mean_x,
    mean_y=kde_stats.mean_y - hist_stats.mean_y,
  )

  alpha = float(meta["alpha"])
  rho = float(meta["rho"])
  bw_x, bw_y = ndkeys_bandwidths(target, rho)
  chi2 = float(meta["chi2"])
  ndf = float(meta["ndf"])
  rchi2 = float(meta.get("reduced_chi2", chi2 / max(ndf, 1)))

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.15)

  x_label = 0.10
  x_int = 0.24
  x_mx = 0.34
  x_my = 0.44
  y_title = 0.95
  y_header = 0.8
  y_hist = 0.55

  latex.SetTextAlign(23)
  latex.DrawLatex(0.28, y_title, "Statistics")
  latex.DrawLatex(x_int, y_header, "integral")
  latex.DrawLatex(x_mx, y_header, "mean x")
  latex.DrawLatex(x_my, y_header, "mean y")

  latex.SetTextAlign(13)
  latex.DrawLatex(x_label, y_hist, "histogram")
  latex.DrawLatex(x_label, y_hist - 0.2, "KDE")
  latex.DrawLatex(x_label, y_hist - 0.4, "difference")

  latex.SetTextAlign(23)
  for row, stats in enumerate((hist_stats, kde_stats, delta)):
    y = y_hist - 0.2 * row
    latex.DrawLatex(x_int, y, f"{stats.integral:.4g}")
    latex.DrawLatex(x_mx, y, f"{stats.mean_x:.4g}")
    latex.DrawLatex(x_my, y, f"{stats.mean_y:.4g}")

  param_lines = [
    f"h_{{x}} = {bw_x:.5g},  h_{{y}} = {bw_y:.5g},  #alpha = {alpha:.5g}",
    f"#chi^{{2}} = {chi2:.4g},  #chi^{{2}}/ndf = {rchi2:.4g},  ndf = {ndf:.0f}",
  ]
  if meta.get("linear_combo", 0):
    mix = float(meta["mix"])
    param_lines.append(
      f"mix = {mix:.4g} (unmirrored),  {1.0 - mix:.4g} (mirrored)"
    )

  y_param = 0.78
  latex.SetTextAlign(12)
  latex.SetTextSize(0.15)
  for line in param_lines:
    latex.DrawLatex(0.58, y_param, line)
    y_param -= 0.20


def plot_overlay(
  target: ROOT.TH2,
  kde_surf: ROOT.TH2,
  meta: dict[str, float | str | int],
  hist_stats: Th2Stats,
  kde_stats: Th2Stats,
  outfile: str,
) -> None:
  data = target.Clone("overlay_data")
  kde_only_surf = kde_surf.Clone("overlay_kde_only_surf")
  kde_only_surf.SetDirectory(0)
  data.SetDirectory(0)

  canvas = ROOT.TCanvas("c_overlay", "2D KDE overlay", 3000, 1640)

  pad_info = ROOT.TPad("pad_info", "", 0.0, 0.0, 1.0, 0.22)
  pad_info.SetFillStyle(0)
  pad_info.Draw()

  pad_left = ROOT.TPad("pad_left", "", 0.0, 0.22, 0.5, 1.0)
  pad_left.Draw()

  pad_right = ROOT.TPad("pad_right", "", 0.5, 0.22, 1.0, 1.0)
  pad_right.Draw()

  pad_left.cd()
  _configure_surf_canvas(pad_left)
  data.SetTitle("Weighted, adaptive, mirrored/unmirrored KDE")
  axis_titles = _style_surf_hist(data, line_color=ROOT.kBlue + 1)
  kde_axes = _style_surf_hist(kde_surf, line_color=ROOT.kRed + 1, line_width=2)
  data.Draw("LEGO")
  kde_surf.SetLineColorAlpha(ROOT.kRed + 1, 0.1)
  kde_surf.Draw("SURF SAME")
  _draw_surf3d_axis_titles(*axis_titles)

  leg = ROOT.TLegend(0.12, 0.82, 0.42, 0.92)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(data, "Data", "l").SetLineWidth(4)
  leg.AddEntry(kde_surf, "#alpha#timesKDE(x,y)", "l").SetLineWidth(4)
  leg.Draw()

  pad_right.cd()
  _configure_surf_canvas(pad_right)
  kde_only_surf.SetTitle("#alpha#timesKDE(x,y)")
  _style_surf_hist(kde_only_surf, line_color=ROOT.kBlue + 1)
  kde_only_surf.Draw("SURF")
  _draw_surf3d_axis_titles(*kde_axes)

  _draw_stats_and_params(pad_info, target, hist_stats, kde_stats, meta)

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def _style_projection_hist(hist: ROOT.TH1) -> None:
  hist.SetLineColor(ROOT.kBlue + 1)
  hist.SetLineWidth(2)


def _draw_projection_panel(
  data: ROOT.TH1,
  model_hist: ROOT.TH1,
  title: str,
) -> None:
  data.SetStats(0)
  data.SetTitle(title)
  data.SetMarkerSize(0.8)
  data.SetMarkerColor(ROOT.kBlack)
  data.SetLineColor(ROOT.kBlack)
  data.SetLineWidth(1)
  model_hist.SetStats(0)
  _style_projection_hist(model_hist)

  ymax = max(data.GetMaximum(), model_hist.GetMaximum())
  ymin = min(0.0, data.GetMinimum())
  span = max(ymax - ymin, 1.0)
  data.GetYaxis().SetRangeUser(ymin - 0.05 * span, ymax + 0.15 * span)
  data.Draw("E1 HIST")
  model_hist.Draw("L SAME")

  leg = ROOT.TLegend(0.62, 0.78, 0.88, 0.90)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.03)
  leg.AddEntry(data, "Data", "lep")
  leg.AddEntry(model_hist, "#alpha#timesKDE", "l")
  leg.Draw()


def _debug_log(hypothesis_id: str, location: str, message: str, data: dict) -> None:
  # #region agent log
  entry = {
    "sessionId": "8940bd",
    "runId": os.environ.get("DEBUG_RUN_ID", "post-fix"),
    "hypothesisId": hypothesis_id,
    "location": location,
    "message": message,
    "data": data,
    "timestamp": int(time.time() * 1000),
  }
  with open(DEBUG_LOG_PATH, "a", encoding="utf-8") as log_file:
    log_file.write(json.dumps(entry) + "\n")
  # #endregion


def plot_projections(
  target: ROOT.TH2,
  kde_model: KdeModel,
  outfile: str,
) -> None:
  if target.GetSumw2N() == 0:
    target.Sumw2()

  proj_opt = "e"
  hx_data = target.ProjectionX(f"{target.GetName()}_px_data", 0, -1, proj_opt)
  hy_data = target.ProjectionY(f"{target.GetName()}_py_data", 0, -1, proj_opt)
  hx_kde = kde_projection_x(kde_model, target, f"{target.GetName()}_px_kde")
  hy_kde = kde_projection_y(kde_model, target, f"{target.GetName()}_py_kde")

  # #region agent log
  x_probe = 0.0
  x_bin = hx_data.GetXaxis().FindBin(x_probe)
  y_bin = hy_data.GetXaxis().FindBin(x_probe)
  _debug_log(
    "H1-H2",
    "plot_2d_kde.py:plot_projections",
    "projection scale check",
    {
      "x_probe": x_probe,
      "data_px": hx_data.GetBinContent(x_bin),
      "kde_px": hx_kde.GetBinContent(hx_kde.GetXaxis().FindBin(x_probe)),
      "data_py": hy_data.GetBinContent(y_bin),
      "kde_py": hy_kde.GetBinContent(hy_kde.GetXaxis().FindBin(x_probe)),
      "sum_data_px": hx_data.Integral(),
      "sum_kde_px": hx_kde.Integral(),
      "sum_data_py": hy_data.Integral(),
      "sum_kde_py": hy_kde.Integral(),
      "n_integrate": KDE_PLOT_BINS,
      "y_scale": target.GetNbinsY() / KDE_PLOT_BINS,
      "x_scale": target.GetNbinsX() / KDE_PLOT_BINS,
    },
  )
  # #endregion

  for h in (hx_data, hy_data, hx_kde, hy_kde):
    h.SetDirectory(0)
  if hx_data.GetSumw2N() == 0:
    hx_data.Sumw2()
  if hy_data.GetSumw2N() == 0:
    hy_data.Sumw2()

  xtitle = target.GetXaxis().GetTitle() or "x [cm]"
  ytitle = target.GetYaxis().GetTitle() or "y [cm]"
  ztitle = target.GetZaxis().GetTitle() or "Entries"
  hx_data.GetXaxis().SetTitle(xtitle)
  hy_data.GetXaxis().SetTitle(ytitle)
  hx_data.GetYaxis().SetTitle(ztitle)
  hy_data.GetYaxis().SetTitle(ztitle)

  canvas = ROOT.TCanvas("c_proj", "2D KDE projections", 2800, 1240)
  canvas.Divide(2, 1)

  canvas.cd(1)
  pad_x = canvas.GetPad(1)
  pad_x.SetGridy()
  pad_x.SetLeftMargin(0.12)
  pad_x.SetBottomMargin(0.12)
  _draw_projection_panel(hx_data, hx_kde, "2D KDE X projection")

  canvas.cd(2)
  pad_y = canvas.GetPad(2)
  pad_y.SetGridy()
  pad_y.SetLeftMargin(0.12)
  pad_y.SetBottomMargin(0.12)
  _draw_projection_panel(hy_data, hy_kde, "2D KDE Y projection")

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def main() -> int:
  show = "--show" in sys.argv
  target, _template, meta, hist_stats, kde_stats = load_fit_objects(FIT_ROOT_FILE)
  rho = float(meta["rho"])
  print_ndkeys_bandwidths(target, rho)
  kde_model = build_kde_model(target, meta)
  kde_fine = kde_plot_hist(kde_model, "kde_fine_plot")
  kde_on_data = evaluate_kde_on_hist_grid(kde_model, target, "kde_on_data_plot")

  os.makedirs(OUTPUT_DIR, exist_ok=True)

  plot_overlay(target, kde_fine, meta, hist_stats, kde_stats, OUTPUT_OVERLAY)
  plot_projections(target, kde_model, OUTPUT_PROJECTIONS)
  plot_ratio(target, kde_on_data, OUTPUT_RATIO)

  if show:
    ROOT.gROOT.SetBatch(False)
    c = ROOT.TCanvas("preview", "preview", 1400, 700)
    c.SetTheta(28)
    c.SetPhi(60)
    c.Divide(2, 1)
    c.cd(1)
    target.SetLineColor(ROOT.kBlue + 1)
    target.SetFillStyle(0)
    kde_surf = kde_fine.Clone("preview_kde_surf")
    kde_surf.SetDirectory(0)
    kde_surf.SetLineColor(ROOT.kRed + 1)
    kde_surf.SetFillStyle(0)
    target.Draw("SURF E")
    kde_surf.Draw("SURF SAME")
    c.cd(2)
    kde_surf.Draw("SURF")
    c.Update()
    input("Press Enter to close...")

  return 0


if __name__ == "__main__":
  sys.exit(main())
