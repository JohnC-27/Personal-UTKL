#!/usr/bin/env python3
"""Plot target data vs weighted adaptive KDE fit from weighted_fixed_kde.root."""

import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

FIT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "weighted_fixed_kde.root"
)
OUTPUT_IMAGE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)),
  "root_files",
  "weighted_fixed_kde_plot.png",
)


def parse_fit_meta(meta: ROOT.TNamed) -> dict[str, float]:
  out = {}
  for part in meta.GetTitle().split(";"):
    key, val = part.split("=", 1)
    out[key] = float(val)
  return out


def load_fit_objects(filepath: str):
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  target = tfile.Get("target_hist")
  template = tfile.Get("kde_template")
  meta = tfile.Get("fit_meta")
  for name, obj in [("target_hist", target), ("kde_template", template), ("fit_meta", meta)]:
    if not obj:
      tfile.Close()
      raise KeyError(f"missing {name!r} in {filepath}")

  target.SetDirectory(0)
  target.SetStats(0)
  template.SetDirectory(0)
  template.SetStats(0)
  meta_title = meta.GetTitle()
  tfile.Close()
  return target, template, parse_fit_meta(meta), meta_title


def _style_target(target: ROOT.TH1) -> None:
  target.SetMarkerSize(0.8)
  target.SetLineWidth(1)


def _style_template(template: ROOT.TH1) -> None:
  template.SetLineColor(ROOT.kBlue + 1)
  template.SetLineWidth(2)
  template.SetFillStyle(0)


def _make_ratio_hist(data: ROOT.TH1, model: ROOT.TH1) -> ROOT.TH1:
  ratio = data.Clone("ratio_data_over_model")
  ratio.SetDirectory(0)
  if ratio.GetSumw2N() == 0:
    ratio.Sumw2()
  ratio.Divide(model)
  ratio.SetTitle("")
  ratio.SetMarkerSize(data.GetMarkerSize())
  ratio.SetMarkerStyle(data.GetMarkerStyle())
  ratio.SetMarkerColor(data.GetMarkerColor())
  ratio.SetLineColor(data.GetLineColor())
  ratio.SetLineWidth(data.GetLineWidth())
  ratio.GetYaxis().SetTitle("Data / Model")
  ratio.GetYaxis().SetTitleOffset(.7)
  ratio.GetYaxis().SetTitleSize(0.07)
  ratio.GetYaxis().SetLabelSize(0.07)
  ratio.GetYaxis().SetNdivisions(505)
  ratio.GetXaxis().SetTitle(data.GetXaxis().GetTitle())
  ratio.GetXaxis().SetLabelSize(0.07)
  ratio.GetXaxis().SetTitleSize(0.07)
  return ratio


def _ratio_y_range(ratio: ROOT.TH1, pad_frac: float = 0.15) -> None:
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


def _draw_unity_line(hist: ROOT.TH1) -> None:
  line = ROOT.TLine(
    hist.GetXaxis().GetXmin(),
    1.0,
    hist.GetXaxis().GetXmax(),
    1.0,
  )
  line.SetLineStyle(2)
  line.SetLineColor(ROOT.kBlack)
  line.Draw()


def plot_fit(
  target: ROOT.TH1,
  template: ROOT.TH1,
  meta: dict[str, float],
  outfile: str,
  show: bool = False,
) -> None:
  alpha = meta["alpha"]
  rho = meta["rho"]
  chi2 = meta["chi2"]
  ndf = meta["ndf"]
  rchi2 = meta.get("reduced_chi2", chi2 / max(ndf, 1))

  _style_target(target)
  _style_template(template)

  canvas = ROOT.TCanvas("c", "Weighted adaptive KDE fit", 1400, 600)
  canvas.Divide(2, 1)

  ratio = _make_ratio_hist(target, template)

  canvas.cd(1)
  pad_left = canvas.GetPad(1)
  pad_left.cd()

  pad_main = ROOT.TPad("pad_main", "", 0.0, 0.32, 1.0, 1.0)
  pad_main.SetBottomMargin(0.02)
  pad_main.SetLeftMargin(0.12)
  pad_main.SetGridy()
  pad_main.Draw()

  pad_ratio = ROOT.TPad("pad_ratio", "", 0.0, 0.0, 1.0, 0.32)
  pad_ratio.SetTopMargin(0.02)
  pad_ratio.SetBottomMargin(0.35)
  pad_ratio.SetLeftMargin(0.12)
  pad_ratio.SetGridy()
  pad_ratio.Draw()

  pad_main.cd()
  target.SetTitle("MM1 X Projection - Weighted Fixed KDE")
  target.GetXaxis().SetLabelSize(0)
  target.GetXaxis().SetTitleSize(0)
  target.Draw("E1 HIST")
  template.Draw("HIST C SAME")

  leg = ROOT.TLegend(0.7, 0.8, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.028)
  leg.AddEntry(target, "Data", "lep")
  leg.AddEntry(template, "#alpha#timesKDE(x)", "l")
  leg.Draw()

  pad_ratio.cd()
  _ratio_y_range(ratio)
  ratio.Draw("E1")
  _draw_unity_line(ratio)

  canvas.cd(2)
  pad_kde = canvas.GetPad(2)
  pad_kde.SetGridy()
  template.SetTitle("Weighted Fixed KDE")
  template.Draw("HIST C")

  leg_kde = ROOT.TLegend(0.4, 0.6, 0.88, 0.88)
  leg_kde.SetBorderSize(0)
  leg_kde.SetFillStyle(0)
  leg_kde.SetTextSize(0.028)
  leg_kde.AddEntry(template, "#alpha#timesKDE(x)", "l")
  leg_kde.Draw()

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.028)
  latex.DrawLatex(0.4, 0.66, f"#rho={rho:.5g}, #alpha={alpha:.5g}")
  latex.DrawLatex(0.4, 0.62, f"#chi^{{2}}={chi2:.3f}, #chi^{{2}}/ndf={rchi2:.3f}")

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")

  if show:
    ROOT.gROOT.SetBatch(False)
    canvas.Draw()
    input("Press Enter to close...")


def main() -> int:
  show = "--show" in sys.argv
  target, template, meta, _ = load_fit_objects(FIT_ROOT_FILE)
  plot_fit(target, template, meta, OUTPUT_IMAGE, show=show)
  return 0


if __name__ == "__main__":
  sys.exit(main())
