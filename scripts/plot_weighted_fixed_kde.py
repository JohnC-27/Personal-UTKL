#!/usr/bin/env python3
"""Plot target data vs weighted adaptive KDE fit from weighted_fixed_kde.root."""

import math
import os
import sys
from dataclasses import dataclass

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

FIT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "weighted_fixed_kde.root"
)
OUTPUT_IMAGE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)),
  "plots",
  "weighted_fixed_kde_plot.png",
)


def parse_fit_meta(meta: ROOT.TNamed) -> dict[str, float]:
  out = {}
  for part in meta.GetTitle().split(";"):
    key, val = part.split("=", 1)
    out[key] = float(val)
  return out


@dataclass
class DistributionStats:
  mean: float
  std: float


def histogram_distribution_stats(hist: ROOT.TH1) -> DistributionStats:
  """Weighted mean and std from bin centers and contents."""
  integral = 0.0
  mean_num = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    weight = hist.GetBinContent(i)
    if weight <= 0:
      continue
    x = hist.GetBinCenter(i)
    integral += weight
    mean_num += weight * x

  if integral <= 0:
    return DistributionStats(0.0, 0.0)

  mean = mean_num / integral
  var_num = 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    weight = hist.GetBinContent(i)
    if weight <= 0:
      continue
    x = hist.GetBinCenter(i)
    diff = x - mean
    var_num += weight * diff * diff

  return DistributionStats(mean=mean, std=math.sqrt(var_num / integral))


def _draw_stats_table(
  pad: ROOT.TPad,
  hist_stats: DistributionStats,
  kde_stats: DistributionStats,
  *,
  unit: str = "cm",
) -> None:
  """Draw mean/std comparison table at the bottom center of a pad."""
  pad.cd()
  delta = DistributionStats(
    mean=kde_stats.mean - hist_stats.mean,
    std=kde_stats.std - hist_stats.std,
  )

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.035)

  x_label = 0.3
  x_mean = 0.50
  x_std = 0.64
  y_title = 0.7
  y_header = 0.61
  y_hist = 0.55
  y_kde = 0.5
  y_delta = 0.45

  latex.SetTextAlign(23)
  latex.DrawLatex(0.50, y_title, f"Statistics ({unit})")

  latex.SetTextAlign(23)
  latex.DrawLatex(x_mean, y_header, "mean")
  latex.DrawLatex(x_std, y_header, "std dev")

  latex.SetTextAlign(13)
  latex.DrawLatex(x_label, y_hist, "histogram")
  latex.DrawLatex(x_label, y_kde, "KDE")
  latex.DrawLatex(x_label, y_delta, "difference")

  latex.SetTextAlign(23)
  latex.DrawLatex(x_mean, y_hist, f"{hist_stats.mean:.4g}")
  latex.DrawLatex(x_std, y_hist, f"{hist_stats.std:.4g}")
  latex.DrawLatex(x_mean, y_kde, f"{kde_stats.mean:.4g}")
  latex.DrawLatex(x_std, y_kde, f"{kde_stats.std:.4g}")
  latex.DrawLatex(x_mean, y_delta, f"{delta.mean:.4g}")
  latex.DrawLatex(x_std, y_delta, f"{delta.std:.4g}")


def load_fit_objects(filepath: str):
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  target = tfile.Get("target_hist")
  template = tfile.Get("kde_template")
  kde_shape = tfile.Get("kde_shape")
  meta = tfile.Get("fit_meta")
  for name, obj in [
    ("target_hist", target),
    ("kde_template", template),
    ("kde_shape", kde_shape),
    ("fit_meta", meta),
  ]:
    if not obj:
      tfile.Close()
      raise KeyError(f"missing {name!r} in {filepath}")

  target.SetDirectory(0)
  target.SetStats(0)
  template.SetDirectory(0)
  template.SetStats(0)
  kde_shape = kde_shape.Clone("kde_shape_plot")
  meta_title = meta.GetTitle()
  tfile.Close()
  return target, template, kde_shape, parse_fit_meta(meta), meta_title


def _style_target(target: ROOT.TH1) -> None:
  target.SetMarkerSize(0.8)
  target.SetLineWidth(1)


def _style_kde_curve(kde_curve: ROOT.TF1) -> None:
  kde_curve.SetLineColor(ROOT.kBlue + 1)
  kde_curve.SetLineWidth(2)


def make_scaled_kde_curve(kde_shape: ROOT.TF1, alpha: float) -> ROOT.TF1:
  """Analytic #alpha#timesKDE(x) for drawing (template stays binned for ratios)."""
  xmin = kde_shape.GetXmin()
  xmax = kde_shape.GetXmax()
  npx = kde_shape.GetNpx()
  if npx <= 0:  # this is just a backup. npx should already be 10000, set in KDE creation script
    npx = 10000

  def scaled(x, _p):
    return alpha * kde_shape.Eval(x[0])

  curve = ROOT.TF1("kde_fit_scaled", scaled, xmin, xmax, 0)
  curve.SetNpx(npx)
  curve._hold_shape = kde_shape
  return curve


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


# pad_frac gives a bit more breathing room to the ratio plot. 0.03 has small scale, 0.15 has larger scale
def _ratio_y_range(ratio: ROOT.TH1, pad_frac: float = 0.03) -> None:
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
  kde_shape: ROOT.TF1,
  meta: dict[str, float],
  outfile: str,
  show: bool = False,
) -> None:
  alpha = meta["alpha"]
  rho = meta["rho"]
  if "bandwidth" in meta:
    bandwidth = meta["bandwidth"]
  elif "h_0" in meta:
    # Older files stored GetFixedWeight() under h_0 (already rho*h_0).
    bandwidth = meta["h_0"]
  else:
    bandwidth = 0.0
  chi2 = meta["chi2"]
  ndf = meta["ndf"]
  rchi2 = meta.get("reduced_chi2", chi2 / max(ndf, 1))

  kde_curve = make_scaled_kde_curve(kde_shape, alpha)
  _style_target(target)
  _style_kde_curve(kde_curve)

  canvas = ROOT.TCanvas("c", "Weighted adaptive KDE fit", 2800, 1200)
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
  target.SetTitle("MM1 X Projection - Weighted, Fixed, Mirrored/Unmirrored KDE")
  target.GetXaxis().SetLabelSize(0)
  target.GetXaxis().SetTitleSize(0)
  target.Draw("E1 HIST")
  kde_curve.Draw("L SAME")

  leg = ROOT.TLegend(0.7, 0.8, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.035)
  leg.AddEntry(target, "Data", "lep")
  leg.AddEntry(kde_curve, "#alpha#timesKDE(x)", "l")
  leg.Draw()

  pad_ratio.cd()
  _ratio_y_range(ratio)
  ratio.Draw("E1")
  _draw_unity_line(ratio)

  hist_stats = histogram_distribution_stats(target)
  kde_stats = histogram_distribution_stats(template)
  pad_left.cd()
  _draw_stats_table(pad_left, hist_stats, kde_stats)

  canvas.cd(2)
  pad_kde = canvas.GetPad(2)
  pad_kde.SetBottomMargin(0.22)
  pad_kde.SetGridy()
  kde_curve.SetTitle("Weighted, Fixed, Mirrored/Unmirrored KDE")
  kde_curve.GetXaxis().SetTitle(target.GetXaxis().GetTitle())
  kde_curve.GetYaxis().SetTitle(target.GetYaxis().GetTitle())
  kde_curve.Draw("L")

  leg_kde = ROOT.TLegend(0.4, 0.6, 0.88, 0.88)
  leg_kde.SetBorderSize(0)
  leg_kde.SetFillStyle(0)
  leg_kde.SetTextSize(0.035)
  leg_kde.AddEntry(kde_curve, "#alpha#timesKDE(x)", "l")
  leg_kde.Draw()

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.035)
  if bandwidth > 0:
    latex.DrawLatex(
      0.3,
      0.66,
      f"Bandwidth={bandwidth:.5g}, #alpha={alpha:.5g}",
    )
  else:
    latex.DrawLatex(0.35, 0.66, f"#rho={rho:.5g}, #alpha={alpha:.5g}")
  latex.DrawLatex(0.35, 0.6, f"#chi^{{2}}={chi2:.3f}, #chi^{{2}}/ndf={rchi2:.3f}")
  if meta.get("linear_combo", 0):
    mix = meta["mix"]
    latex.DrawLatex(
      0.25,
      0.54,
      f"mix={mix:.3g} (unmirrored), {1.0 - mix:.3g} (mirrored)",
    )

  # _draw_stats_table(pad_kde, hist_stats, kde_stats)

  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"Saved {outfile}")

  if show:
    ROOT.gROOT.SetBatch(False)
    canvas.Draw()
    input("Press Enter to close...")


def main() -> int:
  show = "--show" in sys.argv
  target, template, kde_shape, meta, _ = load_fit_objects(FIT_ROOT_FILE)
  plot_fit(target, template, kde_shape, meta, OUTPUT_IMAGE, show=show)
  return 0


if __name__ == "__main__":
  sys.exit(main())
