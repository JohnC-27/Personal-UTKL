#!/usr/bin/env python3
"""
Compare the 2D KDE shape from 2d_kde.root against rebinned NominalxyposMM1 TH2D
histograms in the mz_nominal 2000-bin ROOT files.

The KDE was trained on nominalxyposMM1 in nominal.root (see 2d_kde.py). This script
rebins the high-resolution mz_nominal histograms, profiles alpha against the full
2D histogram, and writes one high-resolution PNG per run with overlay, ratio, and
X/Y projection panels (matching plot_2d_kde.py, without a KDE-only pad).
"""

import array
import json
import math
import os
import sys
import time

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Rebinned grid size (must evenly divide SOURCE_BINS).
N_OUTPUT_BINS = 100

SOURCE_BINS = 2000
HIST_NAME = "NominalxyposMM1"

KDE_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "2d_kde.root"
)
NOMINAL_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "nominal.root"
)
RUN1_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "mz_nominal_2000bin_run1.root"
)
RUN2_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "mz_nominal_2000bin_run2.root"
)
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "plots")

CANVAS_WIDTH = 3200
CANVAS_HEIGHT = 2600
SURF_PLOT_BINS = 100
RATIO_Z_PAD = 1.05
RATIO_Z_MIN_HALF_WIDTH = 0.05


def parse_fit_meta(meta: ROOT.TNamed) -> dict[str, float | str | int]:
  # #region agent log
  _log_path = os.path.join(
    os.path.dirname(__file__), "..", ".cursor", "debug-8dd097.log"
  )
  def _dbg_log(location: str, message: str, data: dict, hypothesis_id: str, run_id: str = "post-fix") -> None:
    with open(_log_path, "a", encoding="utf-8") as log_f:
      log_f.write(
        json.dumps(
          {
            "sessionId": "8dd097",
            "timestamp": int(time.time() * 1000),
            "location": location,
            "message": message,
            "data": data,
            "hypothesisId": hypothesis_id,
            "runId": run_id,
          }
        )
        + "\n"
      )
  # #endregion

  meta_title = meta.GetTitle()
  # #region agent log
  _dbg_log(
    "test_2d_kde_vs_mz_nominal.py:parse_fit_meta:entry",
    "fit_meta title loaded",
    {"meta_title": meta_title, "parts": meta_title.split(";")},
    "A",
  )
  # #endregion

  out: dict[str, float | str | int] = {}
  for part in meta.GetTitle().split(";"):
    key, val = part.split("=", 1)
    # #region agent log
    _dbg_log(
      "test_2d_kde_vs_mz_nominal.py:parse_fit_meta:part",
      "parsing meta part",
      {"key": key, "val": val, "is_pdf": key == "pdf"},
      "B" if key == "pdf" else "C" if key.startswith("ndkeys") else "D",
    )
    # #endregion
    if key == "pdf":
      out[key] = val
      continue
    try:
      num = float(val)
      out[key] = int(num) if num.is_integer() and "." not in val else num
    except ValueError:
      out[key] = val
      # #region agent log
      _dbg_log(
        "test_2d_kde_vs_mz_nominal.py:parse_fit_meta:non_numeric",
        "stored non-numeric meta value as string",
        {"key": key, "val": val},
        "A",
      )
      # #endregion
  # #region agent log
  _dbg_log(
    "test_2d_kde_vs_mz_nominal.py:parse_fit_meta:success",
    "parse_fit_meta completed",
    {"keys": list(out.keys())},
    "A",
  )
  # #endregion
  return out


def load_kde_shape(filepath: str) -> tuple[ROOT.TH2, dict[str, float | str | int]]:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  kde_shape = tfile.Get("kde_shape")
  meta = tfile.Get("fit_meta")
  if not kde_shape or not meta:
    tfile.Close()
    raise KeyError(f"missing kde_shape or fit_meta in {filepath}")

  kde_shape = kde_shape.Clone("kde_shape_test")
  kde_shape.SetDirectory(0)
  meta_dict = parse_fit_meta(meta)
  tfile.Close()
  return kde_shape, meta_dict


def open_histogram2d(filepath: str, hist_name: str) -> ROOT.TH2:
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


def verify_sumw2_errors(hist: ROOT.TH2, label: str) -> bool:
  """Check that GetBinError(ix,iy) equals sqrt(sumw2[bin]) for every bin."""
  if hist.GetSumw2N() == 0:
    print(f"{label}: no Sumw2 array (GetSumw2N() == 0)")
    return False

  sw2 = hist.GetSumw2()
  max_abs_diff = 0.0
  max_rel_diff = 0.0
  n_mismatch = 0
  n_bins = hist.GetNbinsX() * hist.GetNbinsY()

  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      expected = math.sqrt(max(sw2.At(hist.GetBin(ix, iy)), 0.0))
      abs_diff = abs(err - expected)
      if abs_diff > 1e-6 * max(abs(err), abs(expected), 1.0):
        n_mismatch += 1
      max_abs_diff = max(max_abs_diff, abs_diff)
      if expected > 0:
        max_rel_diff = max(max_rel_diff, abs_diff / expected)

  ok = n_mismatch == 0
  status = "PASS" if ok else "FAIL"
  print(
    f"{label}: Sumw2 error check {status} "
    f"({n_bins} bins, mismatches={n_mismatch}, "
    f"max |diff|={max_abs_diff:.3e}, max rel diff={max_rel_diff:.3e})"
  )
  return ok


def rebin_histogram2d(hist: ROOT.TH2, n_output_bins: int, name: str) -> ROOT.TH2:
  if SOURCE_BINS % n_output_bins != 0:
    raise ValueError(
      f"N_OUTPUT_BINS={n_output_bins} must evenly divide SOURCE_BINS={SOURCE_BINS}"
    )
  if hist.GetNbinsX() != SOURCE_BINS or hist.GetNbinsY() != SOURCE_BINS:
    raise ValueError(
      f"expected {SOURCE_BINS}x{SOURCE_BINS} bins before rebinning, got "
      f"{hist.GetNbinsX()}x{hist.GetNbinsY()}"
    )

  factor = SOURCE_BINS // n_output_bins
  rebinned = hist.Rebin2D(factor, factor, name)
  rebinned.SetDirectory(0)
  return rebinned


def optimal_alpha(kde_shape: ROOT.TH2, hist: ROOT.TH2) -> float:
  num, den = 0.0, 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      if err <= 0:
        continue
      observed = hist.GetBinContent(ix, iy)
      shape = kde_shape.GetBinContent(ix, iy)
      w = 1.0 / (err * err)
      num += w * observed * shape
      den += w * shape * shape
  if den <= 0:
    return 1.0
  return num / den


def chi_squared_vs_hist(
  kde_shape: ROOT.TH2,
  hist: ROOT.TH2,
  alpha: float,
) -> tuple[float, int]:
  chi2 = 0.0
  n_used = 0
  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      if err <= 0:
        continue
      observed = hist.GetBinContent(ix, iy)
      expected = alpha * kde_shape.GetBinContent(ix, iy)
      diff = observed - expected
      chi2 += (diff * diff) / (err * err)
      n_used += 1
  return chi2, n_used


def resample_kde_shape(
  kde_shape: ROOT.TH2,
  ref_hist: ROOT.TH2,
  name: str,
) -> ROOT.TH2D:
  """Evaluate the stored KDE shape on the same bin grid as the test histogram."""
  nbx = ref_hist.GetNbinsX()
  nby = ref_hist.GetNbinsY()
  out = ROOT.TH2D(
    name,
    kde_shape.GetTitle(),
    nbx,
    ref_hist.GetXaxis().GetXmin(),
    ref_hist.GetXaxis().GetXmax(),
    nby,
    ref_hist.GetYaxis().GetXmin(),
    ref_hist.GetYaxis().GetXmax(),
  )
  out.SetDirectory(0)

  for ix in range(1, nbx + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    for iy in range(1, nby + 1):
      y = out.GetYaxis().GetBinCenter(iy)
      out.SetBinContent(ix, iy, kde_shape.Interpolate(x, y))

  out._hold_shape = kde_shape
  return out


def kde_template_histogram(
  kde_shape: ROOT.TH2,
  ref_hist: ROOT.TH2,
  alpha: float,
  name: str,
) -> ROOT.TH2D:
  template = kde_shape.Clone(name)
  template.SetDirectory(0)
  template.Scale(alpha)
  xtitle = ref_hist.GetXaxis().GetTitle() or "x [cm]"
  ytitle = ref_hist.GetYaxis().GetTitle() or "y [cm]"
  template.SetTitle(f"#alpha#timesKDE(x,y);{xtitle};{ytitle}")
  return template


def projection_with_errors(
  hist: ROOT.TH2,
  axis: str,
  name: str,
) -> ROOT.TH1:
  if hist.GetSumw2N() == 0:
    hist.Sumw2()
  if axis == "x":
    proj = hist.ProjectionX(name, 0, -1, "e")
  elif axis == "y":
    proj = hist.ProjectionY(name, 0, -1, "e")
  else:
    raise ValueError(f"axis must be 'x' or 'y', got {axis!r}")
  proj.SetDirectory(0)
  if proj.GetSumw2N() == 0:
    proj.Sumw2()
  return proj


def projection_curve(projection: ROOT.TH1, name: str) -> ROOT.TF1:
  xmin = projection.GetXaxis().GetXmin()
  xmax = projection.GetXaxis().GetXmax()

  def interpolated(x, _p):
    return projection.Interpolate(x[0])

  curve = ROOT.TF1(name, interpolated, xmin, xmax, 0)
  curve.SetNpx(2000)
  curve._hold_projection = projection
  return curve


def evaluate_fit(
  label: str,
  data: ROOT.TH2,
  kde_shape: ROOT.TH2,
) -> dict[str, float | ROOT.TH2 | ROOT.TH1 | ROOT.TF1]:
  shape_on_grid = resample_kde_shape(
    kde_shape,
    data,
    f"kde_shape_on_grid_{label}",
  )
  alpha = optimal_alpha(shape_on_grid, data)
  chi2, n_bins = chi_squared_vs_hist(shape_on_grid, data, alpha)
  ndf = max(n_bins - 2, 1)
  reduced = chi2 / ndf
  template = kde_template_histogram(
    shape_on_grid,
    data,
    alpha,
    f"kde_template_{label}",
  )

  data_integral = data.Integral()
  model_integral = template.Integral()
  rel_integral_diff = (
    abs(data_integral - model_integral) / data_integral if data_integral else math.nan
  )

  proj_x_data = projection_with_errors(data, "x", f"data_px_{label}")
  proj_y_data = projection_with_errors(data, "y", f"data_py_{label}")
  proj_x_model = projection_with_errors(template, "x", f"model_px_{label}")
  proj_y_model = projection_with_errors(template, "y", f"model_py_{label}")
  curve_x = projection_curve(proj_x_model, f"curve_px_{label}")
  curve_y = projection_curve(proj_y_model, f"curve_py_{label}")

  print(f"\n{label}:")
  print(f"  bins used in chi2: {n_bins}")
  print(f"  profiled alpha: {alpha:.6g}")
  print(f"  data integral: {data_integral:.6g}")
  print(f"  model integral: {model_integral:.6g}")
  print(f"  |data-model|/data: {rel_integral_diff:.4%}")
  print(f"  chi2 = {chi2:.6g}, ndf = {ndf}, chi2/ndf = {reduced:.6g}")

  return {
    "alpha": alpha,
    "chi2": chi2,
    "ndf": float(ndf),
    "reduced_chi2": reduced,
    "data_integral": data_integral,
    "model_integral": model_integral,
    "rel_integral_diff": rel_integral_diff,
    "template": template,
    "proj_x_data": proj_x_data,
    "proj_y_data": proj_y_data,
    "proj_x_model": proj_x_model,
    "proj_y_model": proj_y_model,
    "curve_x": curve_x,
    "curve_y": curve_y,
  }


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


def _configure_surf_pad(pad: ROOT.TPad) -> None:
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


def interpolate_th2(hist: ROOT.TH2, n_bins: int, name: str) -> ROOT.TH2D:
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()
  ylo = hist.GetYaxis().GetXmin()
  yhi = hist.GetYaxis().GetXmax()
  out = ROOT.TH2D(name, hist.GetTitle(), n_bins, xlo, xhi, n_bins, ylo, yhi)
  out.SetDirectory(0)
  out.SetStats(0)

  for ix in range(1, n_bins + 1):
    x = out.GetXaxis().GetBinCenter(ix)
    for iy in range(1, n_bins + 1):
      y = out.GetYaxis().GetBinCenter(iy)
      out.SetBinContent(ix, iy, hist.Interpolate(x, y))

  return out


def surf_plot_hist(hist: ROOT.TH2, name: str) -> ROOT.TH2D:
  return interpolate_th2(hist, SURF_PLOT_BINS, name)


def _set_diverging_ratio_palette() -> None:
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


def _style_projection_curve(curve: ROOT.TF1) -> None:
  curve.SetLineColor(ROOT.kBlue + 1)
  curve.SetLineWidth(2)


def _draw_projection_on_pad(
  pad: ROOT.TPad,
  data: ROOT.TH1,
  model_curve: ROOT.TF1,
  title: str,
) -> list:
  pad.cd()
  pad.SetGridy()
  pad.SetLeftMargin(0.14)
  pad.SetRightMargin(0.04)
  pad.SetBottomMargin(0.14)
  pad.SetTopMargin(0.08)

  data_plot = data.Clone(f"{data.GetName()}_panel")
  data_plot.SetDirectory(0)
  data_plot.SetStats(0)
  data_plot.SetTitle(title)
  data_plot.SetMarkerSize(0.8)
  data_plot.SetMarkerColor(ROOT.kBlack)
  data_plot.SetLineColor(ROOT.kBlack)
  data_plot.SetLineWidth(1)
  _style_projection_curve(model_curve)

  xmax = data_plot.GetXaxis().GetXmax()
  xmin = data_plot.GetXaxis().GetXmin()
  ymax = max(
    data_plot.GetMaximum(),
    model_curve.GetMaximum(xmin=xmin, xmax=xmax),
  )
  ymin = min(0.0, data_plot.GetMinimum())
  span = max(ymax - ymin, 1.0)
  data_plot.GetYaxis().SetRangeUser(ymin - 0.05 * span, ymax + 0.15 * span)
  data_plot.Draw("E1 HIST")
  model_curve.Draw("L SAME")

  leg = ROOT.TLegend(0.62, 0.78, 0.88, 0.90)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.03)
  leg.AddEntry(data_plot, "Data", "lep")
  leg.AddEntry(model_curve, "#alpha#timesKDE", "l")
  leg.Draw()

  return [data_plot, leg]


def _draw_overlay_on_pad(
  pad: ROOT.TPad,
  data: ROOT.TH2,
  template: ROOT.TH2,
) -> list:
  pad.cd()
  _configure_surf_pad(pad)

  data_plot = data.Clone(f"{data.GetName()}_overlay")
  data_plot.SetDirectory(0)
  data_plot.SetStats(0)
  kde_surf = surf_plot_hist(template, f"{template.GetName()}_overlay_surf")

  data_plot.SetTitle("Data with #alpha#timesKDE(x,y) surface")
  axis_titles = _style_surf_hist(data_plot, line_color=ROOT.kBlue + 1)
  _style_surf_hist(kde_surf, line_color=ROOT.kRed + 1, line_width=2)
  data_plot.Draw("LEGO")
  kde_surf.SetLineColorAlpha(ROOT.kRed + 1, 0.1)
  kde_surf.Draw("SURF SAME")
  _draw_surf3d_axis_titles(*axis_titles)

  leg = ROOT.TLegend(0.12, 0.82, 0.42, 0.92)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(data_plot, "Data", "l")
  leg.AddEntry(kde_surf, "#alpha#timesKDE(x,y)", "l")
  leg.Draw()

  return [data_plot, kde_surf, leg]


def _draw_ratio_on_pad(
  pad: ROOT.TPad,
  data: ROOT.TH2,
  template: ROOT.TH2,
) -> list:
  pad.cd()
  pad.SetRightMargin(0.14)
  pad.SetLeftMargin(0.12)
  pad.SetBottomMargin(0.12)
  pad.SetTopMargin(0.08)

  ratio = kde_over_data_ratio(template, data, f"{data.GetName()}_ratio")
  zmin, zmax = _ratio_color_range(ratio, data)
  xtitle, ytitle, _ztitle = _axis_titles(data)
  ratio.SetTitle(f"KDE / Data;{xtitle};{ytitle};KDE / Data")
  ratio.GetZaxis().SetRangeUser(zmin, zmax)
  ratio.GetZaxis().SetTitle("KDE / Data")
  ratio.GetXaxis().SetTitle(xtitle)
  ratio.GetYaxis().SetTitle(ytitle)

  _set_diverging_ratio_palette()
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

  return [ratio, latex]


def _draw_fit_params(
  pad: ROOT.TPad,
  fit: dict,
  rho: float,
  mix: float,
  linear_combo: bool,
  run_label: str,
  n_output_bins: int,
) -> ROOT.TLatex:
  pad.cd()
  latex = ROOT.TLatex()
  latex.SetNDC(True)
  latex.SetTextFont(42)
  latex.SetTextSize(0.028)
  latex.SetTextAlign(12)
  alpha = fit["alpha"]
  chi2 = fit["chi2"]
  reduced = fit["reduced_chi2"]
  lines = [
    f"{run_label} ({n_output_bins}#times{n_output_bins} bins)",
  ]
  if linear_combo:
    lines.append(f"#rho={rho:.5g}, mix={mix:.5g}, #alpha={alpha:.5g}")
  else:
    lines.append(f"#rho={rho:.5g}, #alpha={alpha:.5g}")
  lines.append(f"#chi^{{2}}={chi2:.2f}, #chi^{{2}}/ndf={reduced:.3f}")
  y = 0.97
  for line in lines:
    latex.DrawLatex(0.02, y, line)
    y -= 0.035
  return latex


def plot_combined_summary(
  data: ROOT.TH2,
  fit: dict,
  run_label: str,
  rho: float,
  mix: float,
  linear_combo: bool,
  n_output_bins: int,
  outfile: str,
) -> None:
  template = fit["template"]
  canvas = ROOT.TCanvas(
    f"c_summary_{run_label}",
    f"2D KDE summary {run_label}",
    CANVAS_WIDTH,
    CANVAS_HEIGHT,
  )

  pad_overlay = ROOT.TPad("pad_overlay", "", 0.0, 0.42, 0.5, 1.0)
  pad_ratio = ROOT.TPad("pad_ratio", "", 0.5, 0.42, 1.0, 1.0)
  pad_x = ROOT.TPad("pad_x", "", 0.0, 0.0, 0.5, 0.42)
  pad_y = ROOT.TPad("pad_y", "", 0.5, 0.0, 1.0, 0.42)
  for pad in (pad_overlay, pad_ratio, pad_x, pad_y):
    pad.Draw()

  keepalive: list = []
  keepalive.extend(_draw_overlay_on_pad(pad_overlay, data, template))
  keepalive.extend(_draw_ratio_on_pad(pad_ratio, data, template))
  keepalive.extend(
    _draw_projection_on_pad(
      pad_x,
      fit["proj_x_data"],
      fit["curve_x"],
      "X projection",
    )
  )
  keepalive.extend(
    _draw_projection_on_pad(
      pad_y,
      fit["proj_y_data"],
      fit["curve_y"],
      "Y projection",
    )
  )
  keepalive.append(
    _draw_fit_params(
      pad_overlay,
      fit,
      rho,
      mix,
      linear_combo,
      run_label,
      n_output_bins,
    )
  )
  canvas._keepalive = keepalive

  canvas.cd()
  canvas.Modified()
  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def output_png_path(run_label: str, n_output_bins: int) -> str:
  return os.path.join(
    OUTPUT_DIR,
    f"2d_kde_vs_mz_nominal_{run_label}_{n_output_bins}x{n_output_bins}bin.png",
  )


def main() -> int:
  if SOURCE_BINS % N_OUTPUT_BINS != 0:
    raise ValueError(
      f"N_OUTPUT_BINS={N_OUTPUT_BINS} must evenly divide SOURCE_BINS={SOURCE_BINS}"
    )

  rebin_factor = SOURCE_BINS // N_OUTPUT_BINS
  print(f"Combining {SOURCE_BINS}x{SOURCE_BINS} bins -> {N_OUTPUT_BINS}x{N_OUTPUT_BINS} "
        f"(factor {rebin_factor})")
  print(f"Loading 2D KDE shape from {KDE_ROOT_FILE}")
  kde_shape, meta = load_kde_shape(KDE_ROOT_FILE)
  rho = meta["rho"]
  mix = meta.get("mix", float("nan"))
  linear_combo = bool(meta.get("linear_combo", 1))
  train_alpha = meta.get("alpha", float("nan"))
  train_chi2 = meta.get("chi2", float("nan"))
  train_ndf = meta.get("ndf", float("nan"))
  print(
    f"Training fit (nominal.root): rho={rho:.6g}, mix={mix:.6g}, alpha={train_alpha:.6g}, "
    f"chi2/ndf={train_chi2 / max(train_ndf, 1):.6g} (ndf={train_ndf:.0f})"
  )

  print(f"\nReference check on nominal.root {HIST_NAME!r} (training target)")
  nominal_ref = open_histogram2d(NOMINAL_ROOT_FILE, "nominalxyposMM1")
  verify_sumw2_errors(nominal_ref, "nominal.root (reference)")
  nominal_fit = evaluate_fit("nominal_ref", nominal_ref, kde_shape)

  datasets = [
    ("run1", RUN1_FILE),
    ("run2", RUN2_FILE),
  ]

  rebinned: dict[str, ROOT.TH2] = {}
  all_sumw2_ok = True

  for label, path in datasets:
    print(f"\nReading {HIST_NAME!r} from {path}")
    raw = open_histogram2d(path, HIST_NAME)
    print(f"  raw: {raw.GetNbinsX()}x{raw.GetNbinsY()} bins")
    all_sumw2_ok &= verify_sumw2_errors(raw, f"{label} ({SOURCE_BINS}x{SOURCE_BINS} bins)")

    rebinned[label] = rebin_histogram2d(
      raw,
      N_OUTPUT_BINS,
      f"{label}_{N_OUTPUT_BINS}x{N_OUTPUT_BINS}",
    )
    all_sumw2_ok &= verify_sumw2_errors(
      rebinned[label],
      f"{label} ({N_OUTPUT_BINS}x{N_OUTPUT_BINS} bins, after rebin)",
    )

  if not all_sumw2_ok:
    print("\nWarning: at least one Sumw2 error check failed.")

  run1_fit = evaluate_fit("run1", rebinned["run1"], kde_shape)
  run2_fit = evaluate_fit("run2", rebinned["run2"], kde_shape)

  print(
    f"\nTraining vs nominal.root: chi2/ndf={nominal_fit['reduced_chi2']:.6g} "
    f"(profiled alpha={nominal_fit['alpha']:.6g})"
  )

  os.makedirs(OUTPUT_DIR, exist_ok=True)
  for run_label, run_fit, run_data in (
    ("run1", run1_fit, rebinned["run1"]),
    ("run2", run2_fit, rebinned["run2"]),
  ):
    plot_combined_summary(
      run_data,
      run_fit,
      run_label,
      rho,
      mix,
      linear_combo,
      N_OUTPUT_BINS,
      output_png_path(run_label, N_OUTPUT_BINS),
    )
  return 0


if __name__ == "__main__":
  sys.exit(main())
