#!/usr/bin/env python3
"""Plot Y projections of NominalxyposMM1 from mz_nominal 2000-bin run1 and run2."""

import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

HIST_NAME = "NominalxyposMM1"
RUN1_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "mz_nominal_2000bin_run1.root"
)
RUN2_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "mz_nominal_2000bin_run2.root"
)
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "plots")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "mm1_y_projection_run1_run2.png")


def load_y_projection(filepath: str, hist_name: str, out_name: str) -> ROOT.TH1:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  hist2d = tfile.Get(hist_name)
  if not hist2d or not hist2d.InheritsFrom("TH2"):
    tfile.Close()
    raise KeyError(f"missing or invalid TH2 {hist_name!r} in {filepath}")

  if hist2d.GetSumw2N() == 0:
    hist2d.Sumw2()

  proj = hist2d.ProjectionY(out_name, 0, -1, "e")
  proj.SetDirectory(0)
  tfile.Close()
  return proj


def style_projection(hist: ROOT.TH1, color: int, title: str) -> None:
  if hist.GetSumw2N() == 0:
    hist.Sumw2()
  hist.SetStats(0)
  hist.SetTitle(title)
  hist.SetLineColor(color)
  hist.SetMarkerColor(color)
  hist.SetLineWidth(1)
  hist.SetMarkerSize(0.8)


def plot_projections(run1: ROOT.TH1, run2: ROOT.TH1, outfile: str) -> None:
  ytitle = run1.GetYaxis().GetTitle() or "Entries"
  xtitle = run1.GetXaxis().GetTitle() or "y [cm]"

  for hist in (run1, run2):
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle(ytitle)

  ymax = max(run1.GetMaximum(), run2.GetMaximum())
  ymin = min(0.0, run1.GetMinimum(), run2.GetMinimum())
  span = max(ymax - ymin, 1.0)
  run1.GetYaxis().SetRangeUser(ymin - 0.05 * span, ymax + 0.15 * span)

  canvas = ROOT.TCanvas("c_mm1_y_proj", "MM1 Y projection", 2000, 1600)
  canvas.SetGridy()
  canvas.SetLeftMargin(0.12)
  canvas.SetRightMargin(0.04)
  canvas.SetBottomMargin(0.12)

  run1.Draw("HIST")
  run2.Draw("HIST SAME")

  leg = ROOT.TLegend(0.66, 0.74, 0.90, 0.90)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.03)
  leg.AddEntry(run1, f"run1 (mean = {run1.GetMean():.4g}cm)", "lep")
  leg.AddEntry(run2, f"run2 (mean = {run2.GetMean():.4g}cm)", "lep")
  leg.Draw()

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def main() -> int:
  show = "--show" in sys.argv

  run1 = load_y_projection(RUN1_FILE, HIST_NAME, "run1_mm1_py")
  run2 = load_y_projection(RUN2_FILE, HIST_NAME, "run2_mm1_py")

  style_projection(run1, ROOT.kBlue + 1, "Y Projection at MM1")
  style_projection(run2, ROOT.kRed + 1, "")

  os.makedirs(OUTPUT_DIR, exist_ok=True)
  plot_projections(run1, run2, OUTPUT_FILE)

  if show:
    ROOT.gROOT.SetBatch(False)
    plot_projections(run1, run2, OUTPUT_FILE)
    input("Press Enter to close...")

  return 0


if __name__ == "__main__":
  sys.exit(main())
