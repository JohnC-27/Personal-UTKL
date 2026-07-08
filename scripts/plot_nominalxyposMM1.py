#!/usr/bin/env python3
"""Plot TH2D nominalxyposMM1 from nominal.root as 2D COLZ and 3D SURF with errors."""

import os
import sys

import ROOT

#from scripts.test_kde_vs_mz_nominal import N_OUTPUT_BINS, SOURCE_BINS

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "Jan2026studies_nominal.root"
)



HIST_NAME = "nominal_xypos_1"
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
OUTPUT_2D = os.path.join(OUTPUT_DIR, "Jan2026_MM1_colz.pdf")
OUTPUT_3D = os.path.join(OUTPUT_DIR, "Jan2026_MM1_surf3d.pdf")
OUTPUT_PROJECTION = os.path.join(OUTPUT_DIR, "Jan2026_MM1_xyproj.pdf")

REBIN = True
N_OUTPUT_BINS_X = 200
N_OUTPUT_BINS_Y = 200

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


def rebin_histogram(
  hist: ROOT.TH2,
  n_output_bins_x: int,
  n_output_bins_y: int,
) -> ROOT.TH2:
  nbx = hist.GetNbinsX()
  nby = hist.GetNbinsY()
  if nbx % n_output_bins_x != 0:
    raise ValueError(
      f"N_OUTPUT_BINS_X={n_output_bins_x} must evenly divide source x bins ({nbx})"
    )
  if nby % n_output_bins_y != 0:
    raise ValueError(
      f"N_OUTPUT_BINS_Y={n_output_bins_y} must evenly divide source y bins ({nby})"
    )

  factor_x = nbx // n_output_bins_x
  factor_y = nby // n_output_bins_y
  rebinned = hist.Rebin2D(factor_x, factor_y, f"{hist.GetName()}_rebinned")
  rebinned.SetDirectory(0)
  return rebinned


def _hist_axis_ranges(hist: ROOT.TH2) -> tuple[float, float]:
  x_range = hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()
  y_range = hist.GetYaxis().GetXmax() - hist.GetYaxis().GetXmin()
  return x_range, y_range


def _style_histogram(hist: ROOT.TH2) -> None:
  hist.SetStats(0)
  hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() or "x [cm]")
  hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() or "y [cm]")
  hist.GetZaxis().SetTitle(hist.GetZaxis().GetTitle() or "Entries")


def _colz_canvas_size(hist: ROOT.TH2, width: int = 900) -> tuple[int, int]:
  x_range, y_range = _hist_axis_ranges(hist)
  if x_range <= 0:
    return width, 800
  height = max(int(width * y_range / x_range), 500)
  return width, height


def plot_colz(hist: ROOT.TH2, outfile: str) -> None:
  _style_histogram(hist)

  width, height = _colz_canvas_size(hist)
  canvas = ROOT.TCanvas("c_colz", "nominalxyposMM1 COLZ", width, height)
  canvas.SetRightMargin(0.14)
  canvas.SetLeftMargin(0.12)
  canvas.SetBottomMargin(0.12)

  hist.Draw("COLZ")
  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def _style_surf3d_axes(hist: ROOT.TH2) -> tuple[str, str, str]:
  """Clear built-in 3D titles; return strings for manual placement below the view."""
  xtitle = hist.GetXaxis().GetTitle() or "x [cm]"
  ytitle = hist.GetYaxis().GetTitle() or "y [cm]"
  ztitle = hist.GetZaxis().GetTitle() or "Entries"
  for axis in (hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()):
    axis.SetTitle("")
    axis.SetLabelSize(0.022)
    axis.SetNdivisions(505)
  return xtitle, ytitle, ztitle


def _draw_surf3d_axis_titles(xtitle: str, ytitle: str, ztitle: str) -> None:
  """Place titles in pad margins so they do not sit on 3D tick labels."""
  latex = ROOT.TLatex()
  latex.SetNDC(True)
  latex.SetTextFont(42)
  latex.SetTextSize(0.032)
  latex.SetTextAlign(22)
  latex.DrawLatex(0.74, 0.05, xtitle)
  latex.DrawLatex(0.26, 0.05, ytitle)
  latex.SetTextAngle(90)
  latex.DrawLatex(0.97, 0.54, ztitle)


def plot_surf3d(hist: ROOT.TH2, outfile: str) -> None:
  h3d = hist.Clone(f"{hist.GetName()}_surf3d")
  h3d.SetDirectory(0)
  _style_histogram(h3d)
  xtitle, ytitle, ztitle = _style_surf3d_axes(h3d)

  h3d.SetLineColor(ROOT.kBlue + 1)
  h3d.SetLineWidth(1)
  h3d.SetFillStyle(0)

  canvas = ROOT.TCanvas("c_surf", "nominalxyposMM1 3D surface", 1200, 950)
  canvas.SetGrid()
  canvas.SetTheta(28)
  canvas.SetPhi(60)
  canvas.SetLeftMargin(0.07)
  canvas.SetRightMargin(0.07)
  canvas.SetBottomMargin(0.05)
  canvas.SetTopMargin(0.1)

  h3d.Draw("SURF E")
  _draw_surf3d_axis_titles(xtitle, ytitle, ztitle)
  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def plot_projection(hist: ROOT.TH2, outfile: str) -> None:
  """Plot both x and y projections as 1D histograms and saves to single file"""
  if hist.GetSumw2N() == 0:
    hist.Sumw2()

  proj_opt = "e"
  hx = hist.ProjectionX(f"{hist.GetName()}_px", 0, -1, proj_opt)
  hy = hist.ProjectionY(f"{hist.GetName()}_py", 0, -1, proj_opt)
  hx.SetDirectory(0)
  hy.SetDirectory(0)

  for h1d in (hx, hy):
    if h1d.GetSumw2N() == 0:
      h1d.Sumw2()

  xtitle = hist.GetXaxis().GetTitle() or "x [cm]"
  ytitle = hist.GetYaxis().GetTitle() or "y [cm]"
  ztitle = hist.GetZaxis().GetTitle() or "Entries"

  for h1d, axis_title, proj_title in (
    (hx, xtitle, "X Projection at MM1"),
    (hy, ytitle, "Y Projection at MM1"),
  ):
    h1d.SetStats(0)
    h1d.SetLineColor(ROOT.kBlue + 1)
    h1d.SetLineWidth(1)
    h1d.SetMarkerSize(0.8)
    h1d.SetMarkerColor(ROOT.kBlue + 1)
    h1d.GetXaxis().SetTitle(axis_title)
    h1d.GetYaxis().SetTitle(ztitle)
    h1d.SetTitle(proj_title)

  canvas = ROOT.TCanvas("c_proj", "nominalxyposMM1 projections", 1400, 600)
  canvas.Divide(2, 1)

  canvas.cd(1)
  pad_x = canvas.GetPad(1)
  pad_x.SetGridy()
  pad_x.SetLeftMargin(0.12)
  pad_x.SetBottomMargin(0.12)
  hx.Draw("E1 HIST")

  canvas.cd(2)
  pad_y = canvas.GetPad(2)
  pad_y.SetGridy()
  pad_y.SetLeftMargin(0.12)
  pad_y.SetBottomMargin(0.12)
  hy.Draw("E1 HIST")

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def main() -> int:
  show = "--show" in sys.argv
  hist = load_histogram(INPUT_ROOT_FILE, HIST_NAME)
  if REBIN:
    hist = rebin_histogram(hist, N_OUTPUT_BINS_X, N_OUTPUT_BINS_Y)

  plot_colz(hist, OUTPUT_2D)
  plot_surf3d(hist, OUTPUT_3D)
  plot_projection(hist, OUTPUT_PROJECTION)

  if show:
    ROOT.gROOT.SetBatch(False)
    c = ROOT.TCanvas("preview", "preview", 1200, 600)
    c.Divide(2, 1)
    c.cd(1)
    hist.Draw("COLZ")
    c.cd(2)
    hist.Draw("SURF E")
    c.Update()
    input("Press Enter to close...")

  return 0


if __name__ == "__main__":
  sys.exit(main())
