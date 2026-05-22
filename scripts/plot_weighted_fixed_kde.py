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

  canvas.cd(1)
  pad_overlay = canvas.GetPad(1)
  pad_overlay.SetGridy()
  target.SetTitle("X projection: data vs weighted adaptive KDE")
  target.Draw("E1 HIST")
  template.Draw("HIST C SAME")

  leg = ROOT.TLegend(0.58, 0.72, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.028)
  leg.AddEntry(target, "Data", "lep")
  leg.AddEntry(template, "#alpha#timesKDE(x)", "l")
  leg.Draw()

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.028)
  latex.DrawLatex(0.58, 0.66, f"#rho={rho:.4g}, #alpha={alpha:.4g}")
  latex.DrawLatex(0.58, 0.62, f"#chi^{{2}}={chi2:.1f}, #chi^{{2}}/ndf={rchi2:.1f}")

  canvas.cd(2)
  pad_kde = canvas.GetPad(2)
  pad_kde.SetGridy()
  template.SetTitle("KDE template")
  template.Draw("HIST C")

  leg_kde = ROOT.TLegend(0.62, 0.82, 0.88, 0.88)
  leg_kde.SetBorderSize(0)
  leg_kde.SetFillStyle(0)
  leg_kde.SetTextSize(0.028)
  leg_kde.AddEntry(template, "#alpha#timesKDE(x)", "l")
  leg_kde.Draw()

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
