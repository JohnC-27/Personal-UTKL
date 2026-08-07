#!/usr/bin/env python3
"""Plot a general TH2 as 2D COLZ, 3D SURF, and XY projections with stats."""

import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# ===============================SCRIPT PARAMS===============================
# Change INPUT, HIST_NAME, OUTPUT FILES, PLOT TITLES

# Entries, mean_x, mean_y, rms_x, rms_y (ROOT computes these from the TH2).
ROOT.gStyle.SetOptStat("neMmRr")
# Means/RMS/stats use in-range bins only (exclude under/overflow).
ROOT.TH1.StatOverflows(False)

INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "ml_beamshift_nX_100um_25bin_corrected.root"
)

HIST_NAME = "ShiftxyposMM1"
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
OUTPUT_2D = os.path.join(OUTPUT_DIR, "ml_beamshift_nX_100um_25bin_corrected_colz.pdf")
OUTPUT_3D = os.path.join(OUTPUT_DIR, "ml_beamshift_nX_100um_25bin_corrected_surf3d.pdf")
OUTPUT_PROJECTION = os.path.join(OUTPUT_DIR, "ml_beamshift_nX_100um_25bin_corrected_xyproj.pdf")

REBIN = False
N_OUTPUT_BINS_X = 200
N_OUTPUT_BINS_Y = 200

SHOW_STATS = True
# TPaveStats corners in NDC (pad fraction): (x1, y1, x2, y2).
# (0,0)=bottom-left of pad, (1,1)=top-right. COLZ default is top-left
# to clear the z-axis palette; projections default to top-right.
STATS_BOX_COLZ = (0.14, 0.68, 0.42, 0.90)
STATS_BOX_PROJ = (0.37, 0.2, 0.67, 0.4)

# "" defaults to histogram name (or title?)
PLOT_TITLE_3D = "-100um x shift MM1"
PLOT_TITLE_PROJ = "-100um x shift MM1" # adds "X/Y Projection"
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


def clear_under_overflow(hist: ROOT.TH1, entries: float | None = None) -> None:
  """Zero under/overflow bins; keep Entries as the TH2 value (not Integral)."""
  n_entries = hist.GetEntries() if entries is None else entries
  if hist.InheritsFrom("TH2"):
    nbx = hist.GetNbinsX()
    nby = hist.GetNbinsY()
    for ix in range(0, nbx + 2):
      for iy in (0, nby + 1):
        hist.SetBinContent(ix, iy, 0.0)
        if hist.GetSumw2N() > 0:
          hist.SetBinError(ix, iy, 0.0)
    for iy in range(0, nby + 2):
      for ix in (0, nbx + 1):
        hist.SetBinContent(ix, iy, 0.0)
        if hist.GetSumw2N() > 0:
          hist.SetBinError(ix, iy, 0.0)
  else:
    nb = hist.GetNbinsX()
    hist.SetBinContent(0, 0.0)
    hist.SetBinContent(nb + 1, 0.0)
    if hist.GetSumw2N() > 0:
      hist.SetBinError(0, 0.0)
      hist.SetBinError(nb + 1, 0.0)
  hist.ResetStats()
  hist.SetEntries(n_entries)


def th2_stats_summary(hist: ROOT.TH2) -> dict[str, float]:
  """Read ROOT-managed TH2 moments; Entries stay the histogram's stored count."""
  n_entries = hist.GetEntries()
  hist.ResetStats()
  hist.SetEntries(n_entries)
  return {
    "entries": float(hist.GetEntries()),
    "integral": float(hist.Integral()),
    "mean_x": float(hist.GetMean(1)),
    "mean_y": float(hist.GetMean(2)),
    "std_x": float(hist.GetStdDev(1)),
    "std_y": float(hist.GetStdDev(2)),
  }


def print_th2_stats(hist: ROOT.TH2, label: str = "") -> None:
  stats = th2_stats_summary(hist)
  prefix = f"{label}: " if label else ""
  print(
    f"{prefix}entries={stats['entries']:.3g}  integral={stats['integral']:.3g}  "
    f"mean=({stats['mean_x']:.3g}, {stats['mean_y']:.3g})  "
    f"std=({stats['std_x']:.3g}, {stats['std_y']:.3g})"
  )


def _hist_axis_ranges(hist: ROOT.TH2) -> tuple[float, float]:
  x_range = hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()
  y_range = hist.GetYaxis().GetXmax() - hist.GetYaxis().GetXmin()
  return x_range, y_range


def _style_histogram(hist: ROOT.TH1, *, show_stats: bool = SHOW_STATS) -> None:
  hist.SetStats(1 if show_stats else 0)
  if hist.InheritsFrom("TH2"):
    hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() or "x [cm]")
    hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() or "y [cm]")
    hist.GetZaxis().SetTitle(hist.GetZaxis().GetTitle() or "Entries")
  else:
    hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() or "x")
    hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() or "Entries")


def _place_stats_box(
  hist: ROOT.TH1,
  box: tuple[float, float, float, float],
) -> None:
  """Move TPaveStats after Draw()+Update; box is (x1, y1, x2, y2) in NDC."""
  if not SHOW_STATS:
    return
  pave = hist.FindObject("stats")
  if not pave:
    return
  x1, y1, x2, y2 = box
  pave.SetX1NDC(x1)
  pave.SetY1NDC(y1)
  pave.SetX2NDC(x2)
  pave.SetY2NDC(y2)
  pave.SetTextSize(0.028)
  
  # OptStat "n" line is TLatex text stored as TNamed title (ROOT statsEditing).
  # Freeze with SetStats(0) so Update/SaveAs does not rebuild and wipe the edit.
  name_line = pave.GetLine(0)
  if name_line:
    name_line.SetTitle("Statistics")
  pave.SetName(f"stats_{hist.GetName()}")
  hist.SetStats(0)
  pave.Draw()


def _colz_canvas_size(hist: ROOT.TH2, width: int = 900) -> tuple[int, int]:
  x_range, y_range = _hist_axis_ranges(hist)
  if x_range <= 0:
    return width, 800
  height = max(int(width * y_range / x_range), 500)
  return width, height


def plot_colz(hist: ROOT.TH2, outfile: str) -> None:
  _style_histogram(hist)

  width, height = _colz_canvas_size(hist)
  canvas = ROOT.TCanvas("c_colz", f"{hist.GetName()} COLZ", width, height)
  canvas.SetRightMargin(0.14)
  canvas.SetLeftMargin(0.12)
  canvas.SetBottomMargin(0.12)
  canvas.SetTopMargin(0.10)

  hist.Draw("COLZ")
  canvas.Update()
  _place_stats_box(hist, STATS_BOX_COLZ)
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
  # SURF view is crowded; keep stats in console / COLZ / projections.
  _style_histogram(h3d, show_stats=False)
  xtitle, ytitle, ztitle = _style_surf3d_axes(h3d)
  if PLOT_TITLE_3D:
    h3d.SetTitle(PLOT_TITLE_3D)

  h3d.SetLineColor(ROOT.kBlue + 1)
  h3d.SetLineWidth(1)
  h3d.SetFillStyle(0)

  canvas = ROOT.TCanvas("LEGO", f"{hist.GetName()} 3D surface", 1200, 950)
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

  # first/last = 1..nbins excludes the other axis under/overflow.
  proj_opt = "e"
  hx = hist.ProjectionX(
    f"{hist.GetName()}_px", 1, hist.GetNbinsY(), proj_opt
  )
  hy = hist.ProjectionY(
    f"{hist.GetName()}_py", 1, hist.GetNbinsX(), proj_opt
  )
  hx.SetDirectory(0)
  hy.SetDirectory(0)

  # Stats-box "Entries" = in-range integral (sum of weights), not GetEntries().
  for h1d in (hx, hy):
    if h1d.GetSumw2N() == 0:
      h1d.Sumw2()
    clear_under_overflow(h1d)
    h1d.SetEntries(h1d.Integral())

  xtitle = hist.GetXaxis().GetTitle() or "x [cm]"
  ytitle = hist.GetYaxis().GetTitle() or "y [cm]"
  ztitle = hist.GetZaxis().GetTitle() or "Entries"

  base = PLOT_TITLE_PROJ or hist.GetName()
  for h1d, axis_title, proj_title in (
    (hx, xtitle, f"{base} X Projection"),
    (hy, ytitle, f"{base} Y Projection"),
  ):
    _style_histogram(h1d)
    h1d.SetLineColor(ROOT.kBlue + 1)
    h1d.SetLineWidth(1)
    h1d.SetMarkerSize(0.8)
    h1d.SetMarkerColor(ROOT.kBlue + 1)
    h1d.GetXaxis().SetTitle(axis_title)
    h1d.GetYaxis().SetTitle(ztitle)
    h1d.SetTitle(proj_title)

  canvas = ROOT.TCanvas("c_proj", f"{hist.GetName()} projections", 1400, 600)
  canvas.Divide(2, 1)

  canvas.cd(1)
  pad_x = canvas.GetPad(1)
  pad_x.SetGridy()
  pad_x.SetLeftMargin(0.12)
  pad_x.SetBottomMargin(0.12)
  pad_x.SetRightMargin(0.08)
  hx.Draw("E1 HIST")
  canvas.Update()
  _place_stats_box(hx, STATS_BOX_PROJ)

  canvas.cd(2)
  pad_y = canvas.GetPad(2)
  pad_y.SetGridy()
  pad_y.SetLeftMargin(0.12)
  pad_y.SetBottomMargin(0.12)
  pad_y.SetRightMargin(0.08)
  hy.Draw("E1 HIST")
  canvas.Update()
  _place_stats_box(hy, STATS_BOX_PROJ)

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def main() -> int:
  show = "--show" in sys.argv
  hist = load_histogram(INPUT_ROOT_FILE, HIST_NAME)
  if REBIN:
    hist = rebin_histogram(hist, N_OUTPUT_BINS_X, N_OUTPUT_BINS_Y)
  n_entries = hist.GetEntries()
  clear_under_overflow(hist, entries=n_entries)

  #print_th2_stats(hist, label=hist.GetName())

  os.makedirs(OUTPUT_DIR, exist_ok=True)
  #plot_colz(hist, OUTPUT_2D)
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
