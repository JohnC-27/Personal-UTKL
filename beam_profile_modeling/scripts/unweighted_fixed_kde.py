#!/usr/bin/env python3
"""Build an unweighted fixed TKDE at rho=0.18 and save a plot."""

import array
import os

import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.SetBatch(True)

INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "nominal.root"
)
TARGET_HIST_NAME = "nominalxyposMM1"
PROJECTION_AXIS = "y"

OUTPUT_PNG = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "unweighted_fixed_kde_plot.png"
)

RHO = .65
TKDE_OPTIONS = "KernelType:Gaussian;Iteration:Fixed;Mirror:noMirror"


def open_histogram(filepath: str, hist_name: str, projection_axis: str) -> ROOT.TH1:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  obj = tfile.Get(hist_name)
  if not obj:
    tfile.Close()
    raise KeyError(f"object {hist_name!r} not found in {filepath}")

  if obj.InheritsFrom("TH2"):
    if projection_axis == "x":
      hist = obj.ProjectionX(f"{hist_name}_px")
    elif projection_axis == "y":
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


def histogram_to_events(hist: ROOT.TH1) -> array.array:
  xs: list[float] = []
  for i in range(1, hist.GetNbinsX() + 1):
    content = hist.GetBinContent(i)
    if content <= 0:
      continue
    n = int(round(content))
    if n <= 0:
      continue
    xs.extend([hist.GetBinCenter(i)] * n)

  if not xs:
    raise ValueError("histogram has no positive bin content")

  return array.array("d", xs)


def scale_factor(kde_func: ROOT.TF1, hist: ROOT.TH1) -> float:
  """Scale KDE density to histogram counts."""
  num, den = 0.0, 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    content = hist.GetBinContent(i)
    if content <= 0:
      continue
    x = hist.GetBinCenter(i)
    num += content
    den += kde_func.Eval(x)
  return num / den if den > 0 else 1.0


def chi_squared_vs_hist(
  kde_func: ROOT.TF1,
  hist: ROOT.TH1,
  alpha: float,
) -> float:
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
  n = 0
  for i in range(1, hist.GetNbinsX() + 1):
    if hist.GetBinError(i) > 0:
      n += 1
  return n


def plot(
  hist: ROOT.TH1,
  kde_func: ROOT.TF1,
  alpha: float,
  rho: float,
  chi2: float,
  ndf: int,
  outfile: str,
) -> None:
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()

  def scaled(x, _p):
    return alpha * kde_func.Eval(x[0])

  curve = ROOT.TF1("kde_scaled", scaled, xlo, xhi, 0)
  curve.SetNpx(kde_func.GetNpx())
  curve.SetLineColor(ROOT.kBlue + 1)
  curve.SetLineWidth(2)
  curve._hold_kde = kde_func

  hist.SetStats(0)
  hist.SetMarkerSize(0.5)
  hist.SetTitle(f"MM1 Y Projection - Unweighted Fixed Unmirrored KDE")
  hist.Draw("E1 HIST")
  curve.Draw("L SAME")

  leg = ROOT.TLegend(0.68, 0.75, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(hist, "Data", "lep")
  leg.AddEntry(curve, "#alpha#timesKDE(x)", "l")
  leg.Draw()

  reduced = chi2 / max(ndf, 1)
  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.028)
  latex.DrawLatex(0.4, 0.70, f"#rho={rho:.2g}, #alpha={alpha:.3g}")
  latex.DrawLatex(0.4, 0.66, f"#chi^{{2}}={chi2:.3f}, #chi^{{2}}/ndf={reduced:.3f}")

  canvas = ROOT.gPad.GetCanvas()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")


def main() -> int:
  print(f"Reading {TARGET_HIST_NAME!r} from {INPUT_ROOT_FILE}")
  hist = open_histogram(INPUT_ROOT_FILE, TARGET_HIST_NAME, PROJECTION_AXIS)

  events = histogram_to_events(hist)
  xlo = hist.GetXaxis().GetXmin()
  xhi = hist.GetXaxis().GetXmax()

  kde = ROOT.TKDE(len(events), events, xlo, xhi, TKDE_OPTIONS, RHO)
  kde_func = kde.GetFunction()
  if kde_func is None:
    raise RuntimeError("TKDE::GetFunction() returned null")
  kde_func.SetNpx(10000)

  alpha = scale_factor(kde_func, hist)
  chi2 = chi_squared_vs_hist(kde_func, hist, alpha)
  ndf = max(count_chi2_bins(hist) - 1, 1)
  reduced = chi2 / ndf
  print(f"Unweighted TKDE: {len(events)} events, rho={RHO}, h_0={kde.GetFixedWeight()}")
  print(f"chi2 = {chi2:.6g}  ndf = {ndf}  chi2/ndf = {reduced:.6g}")

  canvas = ROOT.TCanvas("c", "unweighted fixed KDE", 900, 650)
  canvas.SetGrid()
  canvas.SetLeftMargin(0.12)
  plot(hist, kde_func, alpha, RHO, chi2, ndf, OUTPUT_PNG)
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
