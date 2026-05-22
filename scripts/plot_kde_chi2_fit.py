#!/usr/bin/env python3
"""Plot target data vs KDE chi-squared fit from kde_chi2_fit.root."""

import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

FIT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "kde_chi2_fit.root"
)
OUTPUT_IMAGE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "kde_chi2_fit_plot.png"
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
  template.SetLineWidth(1)
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
  x_scale = meta.get("x_scale", 1.0)
  chi2 = meta["chi2"]
  ndf = meta["ndf"]
  rchi2 = meta.get("reduced_chi2", chi2 / max(ndf, 1))

  _style_target(target)
  _style_template(template)

  canvas = ROOT.TCanvas("c", "KDE chi2 fit", 1400, 600)
  canvas.Divide(2, 1)

  # Left pad: data histogram with KDE overlaid
  canvas.cd(1)
  pad_overlay = canvas.GetPad(1)
  pad_overlay.SetGridy()
  target.SetTitle("X projection KDE")
  target.Draw("E1 HIST")
  template.Draw("HIST C SAME")

  leg = ROOT.TLegend(0.65, 0.85, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.028)
  leg.AddEntry(template, "#alpha#timesKDE(s x)", "l")
  leg.Draw()

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.028)
  latex.DrawLatex(
    0.68, 0.82,
    f"#rho={rho:.3g}, #alpha={alpha:.3g}, #it{{s}}={x_scale:.3g}",
  )
  latex.DrawLatex(
    0.68, 0.78,
    f"#chi^{{2}}={chi2:.1f}, #chi^{{2}}/ndf={rchi2:.1f}",
  )

  # Right pad: KDE template alone
  canvas.cd(2)
  pad_kde = canvas.GetPad(2)
  pad_kde.SetGridy()
  template.SetTitle("KDE template")
  template.Draw("HIST C")

  leg_kde = ROOT.TLegend(0.68, 0.85, 0.88, 0.88)
  leg_kde.SetBorderSize(0)
  leg_kde.SetFillStyle(0)
  leg_kde.SetTextSize(0.028)
  leg_kde.AddEntry(template, "KDE fit (#alpha#timesKDE)", "l")
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
