#!/usr/bin/env python3
"""Ratio plots of two independent TH2s: 2D LEGO plus X/Y projection ratios."""

import array
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# ===============================SCRIPT PARAMS===============================
# Change INPUTS, HIST NAMES, OUTPUT, PLOT TITLE

NUM_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)),
  "root_files",
  "ml_beamshift_nX_100um_100bin_corrected.root",
)
NUM_HIST_NAME = "ShiftxyposMM1"

DEN_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)),
  "root_files",
  "ml_nominal_100bin_corrected.root",
)
DEN_HIST_NAME = "ShiftxyposMM1"

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
OUTPUT_RATIO = os.path.join(
  OUTPUT_DIR, "ml_beamshift_nX_100um_over_nominal_ratio_lego.pdf"
)
OUTPUT_PROJ_RATIO = os.path.join(
  OUTPUT_DIR, "ml_beamshift_nX_100um_over_nominal_proj_ratio.pdf"
)

PLOT_TITLE = "-100um x shift / nominal MM1"
PLOT_TITLE_X_PROJ = f"{PLOT_TITLE} X Projection"
PLOT_TITLE_Y_PROJ = f"{PLOT_TITLE} Y Projection"
Z_TITLE = "Ratio"

RATIO_Z_PAD = 1.05
RATIO_Z_MIN_HALF_WIDTH = 0.05
# ===========================================================================


def load_histogram(filepath: str, hist_name: str) -> ROOT.TH2:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  hist = tfile.Get(hist_name)
  if not hist or not hist.InheritsFrom("TH2"):
    tfile.Close()
    raise KeyError(f"missing or invalid TH2 {hist_name!r} in {filepath}")

  hist.SetDirectory(0)
  tfile.Close()
  return hist


def _assert_compatible(num: ROOT.TH2, den: ROOT.TH2) -> None:
  if num.GetNbinsX() != den.GetNbinsX() or num.GetNbinsY() != den.GetNbinsY():
    raise ValueError(
      f"binning mismatch: num ({num.GetNbinsX()}x{num.GetNbinsY()}) vs "
      f"den ({den.GetNbinsX()}x{den.GetNbinsY()})"
    )
  for axis_name, n_ax, d_ax in (
    ("x", num.GetXaxis(), den.GetXaxis()),
    ("y", num.GetYaxis(), den.GetYaxis()),
  ):
    if abs(n_ax.GetXmin() - d_ax.GetXmin()) > 1e-9 or abs(
      n_ax.GetXmax() - d_ax.GetXmax()
    ) > 1e-9:
      raise ValueError(
        f"{axis_name}-axis range mismatch: "
        f"num [{n_ax.GetXmin()}, {n_ax.GetXmax()}] vs "
        f"den [{d_ax.GetXmin()}, {d_ax.GetXmax()}]"
      )


def make_ratio(num: ROOT.TH1, den: ROOT.TH1, name: str) -> ROOT.TH1:
  """num / den with Sumw2 error propagation (independent histograms)."""
  if num.InheritsFrom("TH2") and den.InheritsFrom("TH2"):
    _assert_compatible(num, den)
  if num.GetSumw2N() == 0:
    num.Sumw2()
  if den.GetSumw2N() == 0:
    den.Sumw2()

  ratio = num.Clone(name)
  ratio.SetDirectory(0)
  if ratio.GetSumw2N() == 0:
    ratio.Sumw2()
  ratio.Divide(den)
  ratio.SetStats(0)
  return ratio


def _project_axis(hist: ROOT.TH2, axis: str, name: str) -> ROOT.TH1:
  if hist.GetSumw2N() == 0:
    hist.Sumw2()
  if axis == "x":
    proj = hist.ProjectionX(name, 1, hist.GetNbinsY(), "e")
  elif axis == "y":
    proj = hist.ProjectionY(name, 1, hist.GetNbinsX(), "e")
  else:
    raise ValueError(f"axis must be 'x' or 'y', got {axis!r}")
  proj.SetDirectory(0)
  if proj.GetSumw2N() == 0:
    proj.Sumw2()
  return proj


def _set_diverging_ratio_palette() -> None:
  """Blue (ratio < 1) -> white (ratio = 1) -> red (ratio > 1)."""
  stops = array.array("d", [0.0, 0.5, 1.0])
  red = array.array("d", [0.0, 1.0, 1.0])
  green = array.array("d", [0.0, 1.0, 0.0])
  blue = array.array("d", [1.0, 1.0, 0.0])
  ROOT.TColor.CreateGradientColorTable(3, stops, red, green, blue, 255)
  ROOT.gStyle.SetNumberContours(255)


def _ratio_color_range(
  ratio: ROOT.TH2,
  den: ROOT.TH2,
  *,
  center: float = 1.0,
  pad: float = RATIO_Z_PAD,
  min_half_width: float = RATIO_Z_MIN_HALF_WIDTH,
) -> tuple[float, float]:
  max_dev = 0.0
  for ix in range(1, den.GetNbinsX() + 1):
    for iy in range(1, den.GetNbinsY() + 1):
      if den.GetBinContent(ix, iy) <= 0:
        continue
      max_dev = max(max_dev, abs(ratio.GetBinContent(ix, iy) - center))

  half = max(max_dev * pad, min_half_width)
  return center - half, center + half


def _ratio_y_range(ratio: ROOT.TH1, pad_frac: float = 0.08) -> None:
  ymax = 1.0
  ymin = 1.0
  for i in range(1, ratio.GetNbinsX() + 1):
    val = ratio.GetBinContent(i)
    err = ratio.GetBinError(i)
    if err <= 0 and val == 0:
      continue
    ymax = max(ymax, val + err)
    ymin = min(ymin, val - err)
  span = max(ymax - ymin, 0.2)
  margin = pad_frac * span
  ratio.GetYaxis().SetRangeUser(ymin - margin, ymax + margin)


def _draw_unity_line(hist: ROOT.TH1) -> ROOT.TLine:
  line = ROOT.TLine(
    hist.GetXaxis().GetXmin(),
    1.0,
    hist.GetXaxis().GetXmax(),
    1.0,
  )
  line.SetLineStyle(2)
  line.SetLineColor(ROOT.kBlack)
  line.Draw()
  return line


def _style_proj_ratio(ratio: ROOT.TH1, axis_title: str, title: str) -> None:
  ratio.SetTitle(title)
  ratio.SetMarkerStyle(20)
  ratio.SetMarkerSize(0.7)
  ratio.SetMarkerColor(ROOT.kBlue + 1)
  ratio.SetLineColor(ROOT.kBlue + 1)
  ratio.SetLineWidth(1)
  ratio.GetXaxis().SetTitle(axis_title)
  ratio.GetXaxis().SetTitleOffset(1.2)
  ratio.GetYaxis().SetTitle(Z_TITLE)
  ratio.GetYaxis().SetTitleOffset(1.2)
  ratio.GetYaxis().SetNdivisions(505)
  _ratio_y_range(ratio)


def plot_ratio_lego(num: ROOT.TH2, den: ROOT.TH2, outfile: str) -> None:
  ratio = make_ratio(num, den, "hist_ratio")
  zmin, zmax = _ratio_color_range(ratio, den)

  xtitle = num.GetXaxis().GetTitle() or "x [cm]"
  ytitle = num.GetYaxis().GetTitle() or "y [cm]"
  ratio.SetTitle("{};{};{};{}".format(PLOT_TITLE, xtitle, ytitle, Z_TITLE))
  ratio.GetZaxis().SetRangeUser(zmin, zmax)
  ratio.GetZaxis().SetTitle(Z_TITLE)
  ratio.GetZaxis().SetTitleOffset(1.5)
  ratio.GetZaxis().SetLabelSize(0.03)

  ratio.GetXaxis().SetTitle(xtitle)
  ratio.GetXaxis().SetTitleOffset(1.8)
  ratio.GetXaxis().SetLabelSize(0.03)

  ratio.GetYaxis().SetTitle(ytitle)
  ratio.GetYaxis().SetTitleOffset(1.8)
  ratio.GetYaxis().SetLabelSize(0.03)

  _set_diverging_ratio_palette()

  canvas = ROOT.TCanvas("c_ratio", "hist ratio", 900, 820)
  canvas.SetRightMargin(0.14)
  canvas.SetLeftMargin(0.12)
  canvas.SetBottomMargin(0.12)
  ratio.Draw("LEGO")

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def plot_projection_ratios(num: ROOT.TH2, den: ROOT.TH2, outfile: str) -> None:
  """Side-by-side X and Y projection ratios (num/den)."""
  _assert_compatible(num, den)
  xtitle = num.GetXaxis().GetTitle() or "x [cm]"
  ytitle = num.GetYaxis().GetTitle() or "y [cm]"

  num_x = _project_axis(num, "x", f"{num.GetName()}_px_num")
  den_x = _project_axis(den, "x", f"{den.GetName()}_px_den")
  num_y = _project_axis(num, "y", f"{num.GetName()}_py_num")
  den_y = _project_axis(den, "y", f"{den.GetName()}_py_den")

  ratio_x = make_ratio(num_x, den_x, "proj_ratio_x")
  ratio_y = make_ratio(num_y, den_y, "proj_ratio_y")
  _style_proj_ratio(ratio_x, xtitle, PLOT_TITLE_X_PROJ)
  _style_proj_ratio(ratio_y, ytitle, PLOT_TITLE_Y_PROJ)

  canvas = ROOT.TCanvas("c_proj_ratio", "projection ratios", 1400, 600)
  canvas.Divide(2, 1)

  canvas.cd(1)
  pad_x = canvas.GetPad(1)
  pad_x.SetGridy()
  pad_x.SetLeftMargin(0.12)
  pad_x.SetBottomMargin(0.12)
  pad_x.SetRightMargin(0.08)
  ratio_x.Draw("HIST")
  pad_x._unity = _draw_unity_line(ratio_x)

  canvas.cd(2)
  pad_y = canvas.GetPad(2)
  pad_y.SetGridy()
  pad_y.SetLeftMargin(0.12)
  pad_y.SetBottomMargin(0.12)
  pad_y.SetRightMargin(0.08)
  ratio_y.Draw("HIST")
  pad_y._unity = _draw_unity_line(ratio_y)

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def main() -> int:
  num = load_histogram(NUM_ROOT_FILE, NUM_HIST_NAME)
  den = load_histogram(DEN_ROOT_FILE, DEN_HIST_NAME)

  os.makedirs(OUTPUT_DIR, exist_ok=True)
  plot_ratio_lego(num, den, OUTPUT_RATIO)
  plot_projection_ratios(num, den, OUTPUT_PROJ_RATIO)
  return 0


if __name__ == "__main__":
  sys.exit(main())
