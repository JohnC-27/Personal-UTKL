#!/usr/bin/env python3
"""
Spline-extrapolate only the top y row (iy 25) of the 2D KDE, using KDE values
in rows iy 1–24 as anchors. All model plots evaluate RooNDKeysPdf directly on
a 200x200 grid (no TH2 bin-center interpolation).
"""

import array
import os
import sys
from dataclasses import dataclass, field

import ROOT

from plot_2d_kde import (
  FIT_ROOT_FILE,
  KdeModel,
  Th2Stats,
  _configure_surf_canvas,
  _draw_projection_panel,
  _draw_surf3d_axis_titles,
  _ratio_color_range,
  _set_diverging_ratio_palette,
  _style_surf_hist,
  histogram_to_weighted_dataset,
  kde_over_data_ratio,
  load_fit_objects,
  make_ndkeys_pdf,
  ndkeys_bandwidths,
)

# Mirror options must match 2d_kde.py (not plot_2d_kde.py) for correct normalization.
NDKEYS_NO_MIRROR = "a"
NDKEYS_MIRROR_BOTH = "am"

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

NOMINAL_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "nominal.root"
)
TARGET_HIST_NAME = "nominalxyposMM1"
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
OUTPUT_PLOT = os.path.join(OUTPUT_DIR, "2d_kde_y_spline_overlay.png")
OUTPUT_PROJECTIONS = os.path.join(OUTPUT_DIR, "2d_kde_y_spline_projections.png")
OUTPUT_RATIO = os.path.join(OUTPUT_DIR, "2d_kde_y_spline_ratio.png")

SPLINE_ROW = 25
N_EVAL_BINS = 200
RATIO_PIXELS_PER_BIN = 4


def build_fit_kde_model(
  target: ROOT.TH2,
  meta: dict[str, float | str | int],
) -> KdeModel:
  ctx = histogram_to_weighted_dataset(target)
  rho = float(meta["rho"])
  alpha = float(meta["alpha"])
  use_linear_combo = bool(meta.get("linear_combo", 0))
  mix = float(meta.get("mix", 1.0))

  if use_linear_combo:
    pdf_unmirrored = make_ndkeys_pdf(
      "y_spline_kde_unmirrored",
      ctx,
      mirror_options=NDKEYS_NO_MIRROR,
      rho=rho,
    )
    pdf_mirrored = make_ndkeys_pdf(
      "y_spline_kde_mirrored",
      ctx,
      mirror_options=NDKEYS_MIRROR_BOTH,
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
    "y_spline_kde_single",
    ctx,
    mirror_options=NDKEYS_MIRROR_BOTH,
    rho=rho,
  )
  return KdeModel(
    ctx=ctx,
    alpha=alpha,
    mix=1.0,
    use_linear_combo=False,
    pdf_single=pdf_single,
  )


@dataclass
class SplineCorrectedKde:
  """RooNDKeysPdf with iy=SPLINE_ROW replaced by per-x cubic spline in y."""

  model: KdeModel
  ref_hist: ROOT.TH2
  norm_scale: float = 1.0
  spline_row: int = SPLINE_ROW
  _y_anchors: list[float] = field(init=False, repr=False)
  _y_spline_min: float = field(init=False, repr=False)

  def __post_init__(self) -> None:
    nby = self.ref_hist.GetNbinsY()
    if self.spline_row < 2 or self.spline_row > nby:
      raise ValueError(f"spline_row={self.spline_row} invalid for {nby} y bins")
    self._y_anchors = [
      self.ref_hist.GetYaxis().GetBinCenter(iy)
      for iy in range(1, self.spline_row)
    ]
    self._y_spline_min = self.ref_hist.GetYaxis().GetBinLowEdge(self.spline_row)

  def _count_scale(self, n_bins_x: int, n_bins_y: int) -> float:
    """Scale PDF point values to per-bin counts on the requested grid."""
    return (self.ref_hist.GetNbinsX() / n_bins_x) * (
      self.ref_hist.GetNbinsY() / n_bins_y
    )

  def kde_at(self, x: float, y: float) -> float:
    return self.model.scaled_at(x, y) * self.norm_scale

  def _spline_at_x(self, x: float) -> ROOT.TSpline3:
    y_arr = array.array("d", self._y_anchors)
    z_arr = array.array(
      "d",
      [self.kde_at(x, y) for y in self._y_anchors],
    )
    return ROOT.TSpline3(
      f"y_spline_x{x:.6g}",
      y_arr,
      z_arr,
      len(y_arr),
      "cs",
      0,
      0,
    )

  def scaled_at(self, x: float, y: float) -> float:
    if y < self._y_spline_min:
      return self.kde_at(x, y)
    return max(self._spline_at_x(x).Eval(y), 0.0)

  def evaluate_th2(
    self,
    *,
    n_bins_x: int,
    n_bins_y: int,
    xlo: float,
    xhi: float,
    ylo: float,
    yhi: float,
    name: str,
    title: str,
  ) -> ROOT.TH2D:
    out = ROOT.TH2D(name, title, n_bins_x, xlo, xhi, n_bins_y, ylo, yhi)
    out.SetDirectory(0)
    out.SetStats(0)
    out._hold_model = self
    count_scale = self._count_scale(n_bins_x, n_bins_y)

    for ix in range(1, n_bins_x + 1):
      x = out.GetXaxis().GetBinCenter(ix)
      spline: ROOT.TSpline3 | None = None
      for iy in range(1, n_bins_y + 1):
        y = out.GetYaxis().GetBinCenter(iy)
        if y >= self._y_spline_min:
          if spline is None:
            spline = self._spline_at_x(x)
          value = max(spline.Eval(y), 0.0)
        else:
          value = self.kde_at(x, y)
        out.SetBinContent(ix, iy, value * count_scale)

    return out

  def evaluate_pure_hist_grid(self, name: str) -> ROOT.TH2D:
    nbx = self.ref_hist.GetNbinsX()
    nby = self.ref_hist.GetNbinsY()
    out = ROOT.TH2D(
      name,
      self.ref_hist.GetTitle(),
      nbx,
      self.ref_hist.GetXaxis().GetXmin(),
      self.ref_hist.GetXaxis().GetXmax(),
      nby,
      self.ref_hist.GetYaxis().GetXmin(),
      self.ref_hist.GetYaxis().GetXmax(),
    )
    out.SetDirectory(0)
    out.SetStats(0)
    out._hold_model = self.model

    for ix in range(1, nbx + 1):
      x = out.GetXaxis().GetBinCenter(ix)
      for iy in range(1, nby + 1):
        y = out.GetYaxis().GetBinCenter(iy)
        out.SetBinContent(ix, iy, self.kde_at(x, y))

    return out

  def evaluate_on_hist_grid(self, name: str) -> ROOT.TH2D:
    corrected = self.evaluate_th2(
      n_bins_x=self.ref_hist.GetNbinsX(),
      n_bins_y=self.ref_hist.GetNbinsY(),
      xlo=self.ref_hist.GetXaxis().GetXmin(),
      xhi=self.ref_hist.GetXaxis().GetXmax(),
      ylo=self.ref_hist.GetYaxis().GetXmin(),
      yhi=self.ref_hist.GetYaxis().GetXmax(),
      name=name,
      title=self.ref_hist.GetTitle(),
    )
    pure = self.evaluate_pure_hist_grid(f"{name}_pure_ref")
    conserve_integral_in_spline_row(pure, corrected, self.spline_row)
    return corrected

  def evaluate_fine_grid(self, name: str) -> ROOT.TH2D:
    corrected = self.evaluate_th2(
      n_bins_x=N_EVAL_BINS,
      n_bins_y=N_EVAL_BINS,
      xlo=self.ref_hist.GetXaxis().GetXmin(),
      xhi=self.ref_hist.GetXaxis().GetXmax(),
      ylo=self.ref_hist.GetYaxis().GetXmin(),
      yhi=self.ref_hist.GetYaxis().GetXmax(),
      name=name,
      title="#alpha#timesKDE(x,y);x [cm];y [cm]",
    )
    pure = self._evaluate_pure_th2(
      n_bins_x=N_EVAL_BINS,
      n_bins_y=N_EVAL_BINS,
      name=f"{name}_pure_ref",
    )
    conserve_integral_above_y(pure, corrected, self._y_spline_min)
    return corrected

  def _evaluate_pure_th2(
    self,
    *,
    n_bins_x: int,
    n_bins_y: int,
    name: str,
  ) -> ROOT.TH2D:
    out = ROOT.TH2D(
      name,
      self.ref_hist.GetTitle(),
      n_bins_x,
      self.ref_hist.GetXaxis().GetXmin(),
      self.ref_hist.GetXaxis().GetXmax(),
      n_bins_y,
      self.ref_hist.GetYaxis().GetXmin(),
      self.ref_hist.GetYaxis().GetXmax(),
    )
    out.SetDirectory(0)
    out.SetStats(0)
    count_scale = self._count_scale(n_bins_x, n_bins_y)

    for ix in range(1, n_bins_x + 1):
      x = out.GetXaxis().GetBinCenter(ix)
      for iy in range(1, n_bins_y + 1):
        y = out.GetYaxis().GetBinCenter(iy)
        out.SetBinContent(ix, iy, self.kde_at(x, y) * count_scale)

    return out


def conserve_integral_in_spline_row(
  reference: ROOT.TH2,
  corrected: ROOT.TH2,
  row: int,
) -> None:
  """Shift mass only in the spline row so the 2D integral matches reference."""
  deficit = reference.Integral() - corrected.Integral()
  if abs(deficit) <= 0:
    return

  nbx = corrected.GetNbinsX()
  row_sum = sum(corrected.GetBinContent(ix, row) for ix in range(1, nbx + 1))
  if row_sum <= 0:
    return

  for ix in range(1, nbx + 1):
    weight = corrected.GetBinContent(ix, row) / row_sum
    corrected.SetBinContent(
      ix,
      row,
      corrected.GetBinContent(ix, row) + deficit * weight,
    )


def conserve_integral_above_y(
  reference: ROOT.TH2,
  corrected: ROOT.TH2,
  y_min: float,
) -> None:
  """Shift mass only at y >= y_min so the 2D integral matches reference."""
  deficit = reference.Integral() - corrected.Integral()
  if abs(deficit) <= 0:
    return

  top_sum = 0.0
  top_bins: list[tuple[int, int]] = []
  for ix in range(1, corrected.GetNbinsX() + 1):
    for iy in range(1, corrected.GetNbinsY() + 1):
      if corrected.GetYaxis().GetBinCenter(iy) < y_min:
        continue
      content = corrected.GetBinContent(ix, iy)
      top_sum += content
      top_bins.append((ix, iy))

  if top_sum <= 0:
    return

  for ix, iy in top_bins:
    weight = corrected.GetBinContent(ix, iy) / top_sum
    corrected.SetBinContent(
      ix,
      iy,
      corrected.GetBinContent(ix, iy) + deficit * weight,
    )


def open_histogram(filepath: str, hist_name: str) -> ROOT.TH2:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  hist = tfile.Get(hist_name)
  if not hist or not hist.InheritsFrom("TH2"):
    tfile.Close()
    raise KeyError(f"missing or invalid TH2 {hist_name!r} in {filepath}")

  hist = hist.Clone(f"{hist_name}_plot")
  hist.SetDirectory(0)
  tfile.Close()
  return hist


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


def _axis_centers(lo: float, hi: float, n_bins: int) -> list[float]:
  width = (hi - lo) / n_bins
  return [lo + (i - 0.5) * width for i in range(1, n_bins + 1)]


def model_projection_x(
  corrected: SplineCorrectedKde,
  ref_hist: ROOT.TH2,
  name: str,
  *,
  n_bins: int = N_EVAL_BINS,
  n_integrate: int = N_EVAL_BINS,
) -> ROOT.TH1D:
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()
  ylo = ref_hist.GetYaxis().GetXmin()
  yhi = ref_hist.GetYaxis().GetXmax()
  out = ROOT.TH1D(name, name, n_bins, xlo, xhi)
  out.SetDirectory(0)
  out.SetStats(0)
  out._hold_model = corrected

  y_values = _axis_centers(ylo, yhi, n_integrate)
  y_scale = ref_hist.GetNbinsY() / n_integrate
  for ix in range(1, n_bins + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    marginal = sum(corrected.scaled_at(x, y) for y in y_values) * y_scale
    out.SetBinContent(ix, marginal)

  return out


def model_projection_y(
  corrected: SplineCorrectedKde,
  ref_hist: ROOT.TH2,
  name: str,
  *,
  n_bins: int = N_EVAL_BINS,
  n_integrate: int = N_EVAL_BINS,
) -> ROOT.TH1D:
  ylo = ref_hist.GetYaxis().GetXmin()
  yhi = ref_hist.GetYaxis().GetXmax()
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()
  out = ROOT.TH1D(name, name, n_bins, ylo, yhi)
  out.SetDirectory(0)
  out.SetStats(0)
  out._hold_model = corrected

  x_values = _axis_centers(xlo, xhi, n_integrate)
  x_scale = ref_hist.GetNbinsX() / n_integrate
  for iy in range(1, n_bins + 1):
    y = out.GetXaxis().GetBinCenter(iy)
    marginal = sum(corrected.scaled_at(x, y) for x in x_values) * x_scale
    out.SetBinContent(iy, marginal)

  return out


def row_chi2_vs_data(
  data: ROOT.TH2,
  model: ROOT.TH2,
  row: int,
) -> tuple[float, int]:
  chi2 = 0.0
  n_bins = 0
  for ix in range(1, data.GetNbinsX() + 1):
    err = data.GetBinError(ix, row)
    if err <= 0:
      continue
    observed = data.GetBinContent(ix, row)
    expected = model.GetBinContent(ix, row)
    diff = observed - expected
    chi2 += (diff * diff) / (err * err)
    n_bins += 1
  return chi2, n_bins


def data_to_fine_grid(hist: ROOT.TH2, n_bins: int, name: str) -> ROOT.TH2D:
  """Split each coarse bin evenly into (n_bins/nbx) x (n_bins/nby) fine bins."""
  nbx = hist.GetNbinsX()
  nby = hist.GetNbinsY()
  if n_bins % nbx != 0 or n_bins % nby != 0:
    raise ValueError(
      f"cannot map {nbx}x{nby} histogram to {n_bins}x{n_bins} evenly"
    )
  fx = n_bins // nbx
  fy = n_bins // nby
  scale = 1.0 / (fx * fy)

  out = ROOT.TH2D(
    name,
    hist.GetTitle(),
    n_bins,
    hist.GetXaxis().GetXmin(),
    hist.GetXaxis().GetXmax(),
    n_bins,
    hist.GetYaxis().GetXmin(),
    hist.GetYaxis().GetXmax(),
  )
  out.SetDirectory(0)
  out.SetStats(0)

  for ix in range(1, nbx + 1):
    for iy in range(1, nby + 1):
      content = hist.GetBinContent(ix, iy)
      err = hist.GetBinError(ix, iy)
      fx0 = (ix - 1) * fx
      fy0 = (iy - 1) * fy
      for subx in range(fx):
        for suby in range(fy):
          out.SetBinContent(fx0 + subx + 1, fy0 + suby + 1, content * scale)
          out.SetBinError(fx0 + subx + 1, fy0 + suby + 1, err * scale)

  return out


def _draw_right_panel(
  pad: ROOT.TPad,
  target: ROOT.TH2,
  hist_stats: Th2Stats,
  kde_stats: Th2Stats,
  meta: dict[str, float | str | int],
  *,
  chi2_before: float,
  chi2_after: float,
  n_spline_bins: int,
) -> None:
  pad.cd()
  pad.SetFillStyle(0)
  pad.SetGrid(0, 0)

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

  x_label = 0.08
  x_int = 0.34
  x_mx = 0.58
  x_my = 0.82

  latex.SetTextSize(0.042)
  latex.SetTextAlign(23)
  latex.DrawLatex(0.50, 0.94, "Statistics")
  latex.DrawLatex(x_int, 0.86, "integral")
  latex.DrawLatex(x_mx, 0.86, "mean x")
  latex.DrawLatex(x_my, 0.86, "mean y")

  latex.SetTextAlign(13)
  latex.SetTextSize(0.038)
  y_hist = 0.76
  row_gap = 0.07
  latex.DrawLatex(x_label, y_hist, "histogram")
  latex.DrawLatex(x_label, y_hist - row_gap, "KDE (spline y)")
  latex.DrawLatex(x_label, y_hist - 2 * row_gap, "difference")

  latex.SetTextAlign(23)
  for row, stats in enumerate((hist_stats, kde_stats, delta)):
    y = y_hist - row_gap * row
    latex.DrawLatex(x_int, y, f"{stats.integral:.4g}")
    latex.DrawLatex(x_mx, y, f"{stats.mean_x:.4g}")
    latex.DrawLatex(x_my, y, f"{stats.mean_y:.4g}")

  latex.SetTextAlign(23)
  latex.SetTextSize(0.042)
  latex.DrawLatex(0.50, 0.48, "Fit parameters")

  param_lines = [
    f"h_{{x}} = {bw_x:.5g},  h_{{y}} = {bw_y:.5g},  #alpha = {alpha:.5g}",
    f"#rho = {rho:.5g}",
    f"#chi^{{2}} = {chi2:.4g},  #chi^{{2}}/ndf = {rchi2:.4g},  ndf = {ndf:.0f}",
  ]
  if meta.get("linear_combo", 0):
    mix = float(meta["mix"])
    param_lines.append(
      f"mix = {mix:.4g} (unmirrored),  {1.0 - mix:.4g} (mirrored)"
    )
  if "pdf" in meta:
    param_lines.append(f"pdf = {meta['pdf']}")

  param_lines.extend([
    f"y-spline: iy {SPLINE_ROW} only (from iy {SPLINE_ROW - 1})",
    f"eval grid: {N_EVAL_BINS}#times{N_EVAL_BINS} PDF points",
    f"iy {SPLINE_ROW} #chi^{{2}}: {chi2_before:.4g} #rightarrow {chi2_after:.4g} ({n_spline_bins} bins)",
  ])

  latex.SetTextAlign(12)
  latex.SetTextSize(0.034)
  y_param = 0.40
  for line in param_lines:
    latex.DrawLatex(0.08, y_param, line)
    y_param -= 0.055


def plot_overlay_with_stats(
  target: ROOT.TH2,
  kde_surf: ROOT.TH2,
  meta: dict[str, float | str | int],
  hist_stats: Th2Stats,
  kde_stats: Th2Stats,
  spline_info: dict[str, float | int],
  outfile: str,
) -> None:
  data = target.Clone("overlay_data")
  kde_plot = kde_surf.Clone("overlay_kde_surf")
  data.SetDirectory(0)
  kde_plot.SetDirectory(0)
  data.SetStats(0)

  canvas = ROOT.TCanvas("c_y_spline_overlay", "2D KDE y-spline overlay", 2400, 1200)

  pad_left = ROOT.TPad("pad_left", "", 0.0, 0.0, 0.62, 1.0)
  pad_left.Draw()

  pad_right = ROOT.TPad("pad_right", "", 0.62, 0.0, 1.0, 1.0)
  pad_right.Draw()

  pad_left.cd()
  _configure_surf_canvas(pad_left)
  data.SetTitle("Data with spline-corrected #alpha#timesKDE(x,y)")
  axis_titles = _style_surf_hist(data, line_color=ROOT.kBlue + 1)
  _style_surf_hist(kde_plot, line_color=ROOT.kRed + 1, line_width=2)
  data.Draw("LEGO")
  kde_plot.SetLineColorAlpha(ROOT.kRed + 1, 0.12)
  kde_plot.Draw("SURF SAME")
  _draw_surf3d_axis_titles(*axis_titles)

  leg = ROOT.TLegend(0.12, 0.82, 0.50, 0.92)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(data, TARGET_HIST_NAME, "l")
  leg.AddEntry(kde_plot, "KDE + y-spline (iy 25)", "l")
  leg.Draw()

  _draw_right_panel(
    pad_right,
    target,
    hist_stats,
    kde_stats,
    meta,
    chi2_before=float(spline_info["chi2_before"]),
    chi2_after=float(spline_info["chi2_after"]),
    n_spline_bins=int(spline_info["n_spline_bins"]),
  )

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def plot_projections(
  target: ROOT.TH2,
  corrected: SplineCorrectedKde,
  outfile: str,
) -> None:
  if target.GetSumw2N() == 0:
    target.Sumw2()

  proj_opt = "e"
  hx_data = target.ProjectionX(f"{target.GetName()}_px_data", 0, -1, proj_opt)
  hy_data = target.ProjectionY(f"{target.GetName()}_py_data", 0, -1, proj_opt)
  hx_kde = model_projection_x(corrected, target, f"{target.GetName()}_px_kde")
  hy_kde = model_projection_y(corrected, target, f"{target.GetName()}_py_kde")

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

  canvas = ROOT.TCanvas("c_y_spline_proj", "2D KDE y-spline projections", 2800, 1240)
  canvas.Divide(2, 1)

  canvas.cd(1)
  pad_x = canvas.GetPad(1)
  pad_x.SetGridy()
  pad_x.SetLeftMargin(0.12)
  pad_x.SetBottomMargin(0.12)
  _draw_projection_panel(hx_data, hx_kde, "X projection")

  canvas.cd(2)
  pad_y = canvas.GetPad(2)
  pad_y.SetGridy()
  pad_y.SetLeftMargin(0.12)
  pad_y.SetBottomMargin(0.12)
  _draw_projection_panel(hy_data, hy_kde, "Y projection")

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def _ratio_value_to_color(value: float, zmin: float, zmax: float) -> int:
  if zmax <= zmin:
    t = 0.5
  else:
    t = (value - zmin) / (zmax - zmin)
  t = max(0.0, min(1.0, t))
  return ROOT.gStyle.GetColorPalette(int(round(t * 255)))


def plot_ratio_pixels(
  target: ROOT.TH2,
  template: ROOT.TH2,
  outfile: str,
  *,
  pixels_per_bin: int = RATIO_PIXELS_PER_BIN,
) -> None:
  del pixels_per_bin  # kept for API compatibility; bins are drawn at native width

  ratio = kde_over_data_ratio(template, target, "kde_over_data_ratio_pixels")
  zmin, zmax = _ratio_color_range(ratio, target)

  xtitle = target.GetXaxis().GetTitle() or "x [cm]"
  ytitle = target.GetYaxis().GetTitle() or "y [cm]"
  ratio.SetTitle("KDE / Data;{};{};KDE / Data".format(xtitle, ytitle))
  ratio.GetZaxis().SetRangeUser(zmin, zmax)
  ratio.GetZaxis().SetTitle("KDE / Data")
  ratio.GetXaxis().SetTitle(xtitle)
  ratio.GetYaxis().SetTitle(ytitle)

  _set_diverging_ratio_palette()
  ROOT.gStyle.SetCanvasPreferGL(False)

  nbx = ratio.GetNbinsX()
  nby = ratio.GetNbinsY()
  margin_left = 110
  margin_right = 130
  margin_bottom = 110
  margin_top = 70
  canvas_w = 900
  canvas_h = 820
  palette_frac = margin_right / canvas_w
  plot_right = 1.0 - palette_frac

  canvas = ROOT.TCanvas("c_ratio_pixels", "KDE / Data ratio", canvas_w, canvas_h)
  canvas.SetGrid(0, 0)

  plot_pad = ROOT.TPad(
    "plot_pad",
    "",
    margin_left / canvas_w,
    margin_bottom / canvas_h,
    plot_right,
    1.0 - margin_top / canvas_h,
  )
  plot_pad.SetGrid(0, 0)
  plot_pad.Draw()

  palette_pad = ROOT.TPad(
    "palette_pad",
    "",
    plot_right,
    margin_bottom / canvas_h,
    1.0,
    1.0 - margin_top / canvas_h,
  )
  palette_pad.SetGrid(0, 0)
  palette_pad.Draw()

  plot_pad.cd()
  frame = ratio.Clone("ratio_axis_frame")
  frame.SetDirectory(0)
  frame.Reset()
  frame.Draw("AXIS")

  boxes: list[ROOT.TBox] = []
  for ix in range(1, nbx + 1):
    xlo = ratio.GetXaxis().GetBinLowEdge(ix)
    xhi = ratio.GetXaxis().GetBinUpEdge(ix)
    for iy in range(1, nby + 1):
      ylo = ratio.GetYaxis().GetBinLowEdge(iy)
      yhi = ratio.GetYaxis().GetBinUpEdge(iy)
      if target.GetBinContent(ix, iy) <= 0:
        fill_color = ROOT.kWhite
      else:
        fill_color = _ratio_value_to_color(ratio.GetBinContent(ix, iy), zmin, zmax)

      box = ROOT.TBox(xlo, ylo, xhi, yhi)
      box.SetFillColor(fill_color)
      box.SetFillStyle(1001)
      box.SetLineWidth(0)
      box.SetLineColor(fill_color)
      box.Draw("SAME")
      boxes.append(box)

  frame._hold_boxes = boxes

  palette_pad.cd()
  palette_hist = ratio.Clone("ratio_palette_hist")
  palette_hist.SetDirectory(0)
  palette_hist.GetXaxis().SetLabelSize(0)
  palette_hist.GetXaxis().SetTitleSize(0)
  palette_hist.GetYaxis().SetLabelSize(0)
  palette_hist.GetYaxis().SetTitleSize(0)
  palette_hist.GetYaxis().SetTickLength(0)
  palette_hist.Draw("COLZ")

  canvas.cd()
  latex = ROOT.TLatex()
  latex.SetNDC(True)
  latex.SetTextFont(42)
  latex.SetTextSize(0.03)
  latex.SetTextAlign(22)
  latex.DrawLatex(
    0.50,
    0.96,
    "blue: KDE < data   white: KDE = data   red: KDE > data",
  )

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def normalization_scale(
  kde_template: ROOT.TH2,
  model: KdeModel,
  ref_hist: ROOT.TH2,
) -> float:
  """Match live PDF integral on the histogram grid to the fitted template."""
  probe = SplineCorrectedKde(model=model, ref_hist=ref_hist, norm_scale=1.0)
  pure = probe.evaluate_pure_hist_grid("norm_probe_pure")
  template_integral = kde_template.Integral()
  pure_integral = pure.Integral()
  if pure_integral <= 0:
    return 1.0
  return template_integral / pure_integral


def main() -> int:
  target_from_fit, kde_template, meta, _hist_stats_orig, _kde_stats_orig = load_fit_objects(
    FIT_ROOT_FILE
  )
  target = open_histogram(NOMINAL_ROOT_FILE, TARGET_HIST_NAME)
  kde_model = build_fit_kde_model(target_from_fit, meta)
  norm_scale = normalization_scale(kde_template, kde_model, target)
  corrected = SplineCorrectedKde(
    model=kde_model,
    ref_hist=target,
    norm_scale=norm_scale,
  )

  if target.GetSumw2N() == 0:
    target.Sumw2()

  kde_pure = corrected.evaluate_pure_hist_grid("kde_pure_hist_grid")
  kde_on_grid = corrected.evaluate_on_hist_grid("kde_spline_hist_grid")

  chi2_before, n_spline_bins = row_chi2_vs_data(target, kde_pure, SPLINE_ROW)
  chi2_after, _ = row_chi2_vs_data(target, kde_on_grid, SPLINE_ROW)

  hist_stats = th2_distribution_stats(target)
  kde_stats = th2_distribution_stats(kde_on_grid)
  kde_surf = corrected.evaluate_fine_grid("kde_y_spline_surf")

  data_fine = data_to_fine_grid(target, N_EVAL_BINS, "data_fine_ratio")

  os.makedirs(OUTPUT_DIR, exist_ok=True)
  plot_overlay_with_stats(
    target,
    kde_surf,
    meta,
    hist_stats,
    kde_stats,
    spline_info={
      "chi2_before": chi2_before,
      "chi2_after": chi2_after,
      "n_spline_bins": n_spline_bins,
    },
    outfile=OUTPUT_PLOT,
  )
  plot_projections(target, corrected, OUTPUT_PROJECTIONS)
  plot_ratio_pixels(data_fine, kde_surf, OUTPUT_RATIO)

  print(
    f"iy {SPLINE_ROW} chi2: {chi2_before:.4g} -> {chi2_after:.4g} "
    f"over {n_spline_bins} bins  (model eval {N_EVAL_BINS}x{N_EVAL_BINS})"
  )
  return 0


if __name__ == "__main__":
  sys.exit(main())
