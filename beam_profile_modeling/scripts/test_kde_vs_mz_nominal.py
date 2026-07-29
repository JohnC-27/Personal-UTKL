#!/usr/bin/env python3
"""
Compare the fixed KDE from weighted_fixed_kde.root against rebinned X projections
of NominalxyposMM1 in the mz_nominal 2000-bin ROOT files.
"""

import math
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Target number of bins after combining the 2000-bin projections (must divide 2000).
N_OUTPUT_BINS = 100

SOURCE_BINS = 2000
HIST_NAME = "NominalxyposMM1"
PROJECTION_AXIS = "x"

KDE_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "weighted_fixed_kde.root"
)
RUN1_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "mz_nominal_2000bin_run1.root"
)
RUN2_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", "mz_nominal_2000bin_run2.root"
)
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "plots")


def parse_fit_meta(meta: ROOT.TNamed) -> dict[str, float]:
  out: dict[str, float] = {}
  for part in meta.GetTitle().split(";"):
    key, val = part.split("=", 1)
    out[key] = float(val)
  return out


def load_kde_fit(filepath: str) -> tuple[ROOT.TF1, dict[str, float]]:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  kde_shape = tfile.Get("kde_shape")
  meta = tfile.Get("fit_meta")
  if not kde_shape or not meta:
    tfile.Close()
    raise KeyError(f"missing kde_shape or fit_meta in {filepath}")

  kde_shape = kde_shape.Clone("kde_shape_test")
  meta_dict = parse_fit_meta(meta)
  tfile.Close()
  return kde_shape, meta_dict


def open_projected_histogram(
  filepath: str,
  hist_name: str,
  projection_axis: str = "x",
) -> ROOT.TH1:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  obj = tfile.Get(hist_name)
  if not obj:
    tfile.Close()
    raise KeyError(f"object {hist_name!r} not found in {filepath}")

  if obj.InheritsFrom("TH2"):
    axis = projection_axis.lower()
    suffix = "px" if axis == "x" else "py"
    if axis not in ("x", "y"):
      tfile.Close()
      raise ValueError(f"projection_axis must be 'x' or 'y', got {projection_axis!r}")
    hist = obj.ProjectionX(f"{hist_name}_{suffix}", 0, -1, "e")
    if axis == "y":
      hist = obj.ProjectionY(f"{hist_name}_{suffix}", 0, -1, "e")
  elif obj.InheritsFrom("TH1"):
    hist = obj
  else:
    tfile.Close()
    raise TypeError(f"{hist_name!r} is not TH1/TH2")

  hist.SetDirectory(0)
  tfile.Close()
  return hist


def verify_sumw2_errors(hist: ROOT.TH1, label: str) -> bool:
  """Check that GetBinError(i) equals sqrt(sumw2[i]) for every bin."""
  if hist.GetSumw2N() == 0:
    print(f"{label}: no Sumw2 array (GetSumw2N() == 0)")
    return False

  sw2 = hist.GetSumw2()
  max_abs_diff = 0.0
  max_rel_diff = 0.0
  n_mismatch = 0

  for i in range(1, hist.GetNbinsX() + 1):
    err = hist.GetBinError(i)
    expected = math.sqrt(max(sw2.At(i), 0.0))
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
    f"({hist.GetNbinsX()} bins, mismatches={n_mismatch}, "
    f"max |diff|={max_abs_diff:.3e}, max rel diff={max_rel_diff:.3e})"
  )
  return ok


def rebin_histogram(hist: ROOT.TH1, n_output_bins: int, name: str) -> ROOT.TH1:
  if SOURCE_BINS % n_output_bins != 0:
    raise ValueError(
      f"N_OUTPUT_BINS={n_output_bins} must evenly divide SOURCE_BINS={SOURCE_BINS}"
    )
  if hist.GetNbinsX() != SOURCE_BINS:
    raise ValueError(
      f"expected {SOURCE_BINS} bins before rebinning, got {hist.GetNbinsX()}"
    )

  factor = SOURCE_BINS // n_output_bins
  rebinned = hist.Rebin(factor, name)
  rebinned.SetDirectory(0)
  return rebinned


def optimal_alpha(kde_func: ROOT.TF1, hist: ROOT.TH1) -> float:
  num, den = 0.0, 0.0
  for i in range(1, hist.GetNbinsX() + 1):
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    observed = hist.GetBinContent(i)
    shape = kde_func.Eval(hist.GetBinCenter(i))
    w = 1.0 / (err * err)
    num += w * observed * shape
    den += w * shape * shape
  if den <= 0:
    return 1.0
  return num / den


def chi_squared_vs_hist(
  kde_func: ROOT.TF1,
  hist: ROOT.TH1,
  alpha: float,
) -> tuple[float, int]:
  chi2 = 0.0
  n_used = 0
  for i in range(1, hist.GetNbinsX() + 1):
    err = hist.GetBinError(i)
    if err <= 0:
      continue
    observed = hist.GetBinContent(i)
    expected = alpha * kde_func.Eval(hist.GetBinCenter(i))
    diff = observed - expected
    chi2 += (diff * diff) / (err * err)
    n_used += 1
  return chi2, n_used


def kde_template_histogram(
  kde_func: ROOT.TF1,
  ref_hist: ROOT.TH1,
  alpha: float,
  name: str,
) -> ROOT.TH1D:
  nb = ref_hist.GetNbinsX()
  xlo = ref_hist.GetXaxis().GetXmin()
  xhi = ref_hist.GetXaxis().GetXmax()
  xtitle = ref_hist.GetXaxis().GetTitle() or "x [cm]"
  out = ROOT.TH1D(name, f"#alpha#timesKDE(x);{xtitle};Entries", nb, xlo, xhi)
  out.SetDirectory(0)
  for i in range(1, nb + 1):
    x = ref_hist.GetBinCenter(i)
    out.SetBinContent(i, alpha * kde_func.Eval(x))
  return out


def make_scaled_kde_curve(kde_shape: ROOT.TF1, alpha: float, name: str) -> ROOT.TF1:
  xmin = kde_shape.GetXmin()
  xmax = kde_shape.GetXmax()
  npx = kde_shape.GetNpx()
  if npx <= 0:
    npx = 10000

  def scaled(x, _p):
    return alpha * kde_shape.Eval(x[0])

  curve = ROOT.TF1(name, scaled, xmin, xmax, 0)
  curve.SetNpx(npx)
  curve._hold_shape = kde_shape
  return curve


def evaluate_fit(
  label: str,
  data: ROOT.TH1,
  kde_shape: ROOT.TF1,
) -> dict[str, float | ROOT.TH1 | ROOT.TF1]:
  alpha = optimal_alpha(kde_shape, data)
  chi2, n_bins = chi_squared_vs_hist(kde_shape, data, alpha)
  ndf = max(n_bins - 1, 1)
  reduced = chi2 / ndf
  template = kde_template_histogram(kde_shape, data, alpha, f"kde_template_{label}")
  curve = make_scaled_kde_curve(kde_shape, alpha, f"kde_curve_{label}")

  data_integral = data.Integral()
  model_integral = template.Integral()
  rel_integral_diff = (
    abs(data_integral - model_integral) / data_integral if data_integral else math.nan
  )

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
    "curve": curve,
  }


def _make_ratio_hist(data: ROOT.TH1, model: ROOT.TH1, name: str) -> ROOT.TH1:
  ratio = data.Clone(name)
  ratio.SetDirectory(0)
  if ratio.GetSumw2N() == 0:
    ratio.Sumw2()
  ratio.Divide(model)
  ratio.SetTitle("")
  ratio.GetYaxis().SetTitle("Data / Model")
  ratio.GetYaxis().SetTitleOffset(0.7)
  ratio.GetYaxis().SetTitleSize(0.07)
  ratio.GetYaxis().SetLabelSize(0.07)
  ratio.GetYaxis().SetNdivisions(505)
  ratio.GetXaxis().SetTitle(data.GetXaxis().GetTitle())
  ratio.GetXaxis().SetLabelSize(0.07)
  ratio.GetXaxis().SetTitleSize(0.07)
  return ratio


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


def _draw_run_panel(
  pad_index: int,
  canvas: ROOT.TCanvas,
  data: ROOT.TH1,
  fit: dict,
  panel_title: str,
  rho: float,
) -> list:
  canvas.cd(pad_index)
  pad_outer = canvas.GetPad(pad_index)
  pad_outer.cd()

  pad_main = ROOT.TPad(f"pad_main_{pad_index}", "", 0.0, 0.32, 1.0, 1.0)
  pad_main.SetBottomMargin(0.02)
  pad_main.SetLeftMargin(0.14)
  pad_main.SetRightMargin(0.04)
  pad_main.SetGridy()
  pad_main.Draw()

  pad_ratio = ROOT.TPad(f"pad_ratio_{pad_index}", "", 0.0, 0.0, 1.0, 0.32)
  pad_ratio.SetTopMargin(0.02)
  pad_ratio.SetBottomMargin(0.35)
  pad_ratio.SetLeftMargin(0.14)
  pad_ratio.SetRightMargin(0.04)
  pad_ratio.SetGridy()
  pad_ratio.Draw()

  curve = fit["curve"]
  template = fit["template"]
  alpha = fit["alpha"]
  chi2 = fit["chi2"]
  reduced = fit["reduced_chi2"]

  data_plot = data.Clone(f"data_panel_{pad_index}")
  data_plot.SetDirectory(0)
  data_plot.SetStats(0)
  data_plot.SetMarkerSize(0.7)
  data_plot.SetLineWidth(1)

  curve.SetLineColor(ROOT.kBlue + 1)
  curve.SetLineWidth(2)

  ratio = _make_ratio_hist(data_plot, template, f"ratio_{pad_index}")

  pad_main.cd()
  data_plot.SetTitle(panel_title)
  data_plot.GetXaxis().SetLabelSize(0)
  data_plot.GetXaxis().SetTitleSize(0)
  data_plot.Draw("E1 HIST")
  curve.Draw("L SAME")

  leg = ROOT.TLegend(0.75, 0.72, 0.92, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.03)
  leg.AddEntry(data_plot, "Data", "lep")
  leg.AddEntry(curve, "#alpha#timesKDE(x)", "l")
  leg.Draw()

  pad_main.cd()
  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextFont(42)
  latex.SetTextSize(0.03)
  latex.DrawLatex(0.16, 0.82, f"#rho={rho:.5g}, #alpha={alpha:.5g}")
  latex.DrawLatex(
    0.16,
    0.76,
    f"#chi^{{2}}={chi2:.2f}, #chi^{{2}}/ndf={reduced:.3f}",
  )
  pad_main.Modified()
  pad_main.Update()

  pad_ratio.cd()
  _ratio_y_range(ratio)
  ratio.Draw("E1")
  unity = ROOT.TLine(
    ratio.GetXaxis().GetXmin(),
    1.0,
    ratio.GetXaxis().GetXmax(),
    1.0,
  )
  unity.SetLineStyle(2)
  unity.SetLineColor(ROOT.kBlack)
  unity.Draw()
  pad_ratio.Modified()
  pad_ratio.Update()

  pad_outer.cd()
  pad_outer.Modified()
  pad_outer.Update()

  return [
    pad_outer,
    pad_main,
    pad_ratio,
    data_plot,
    curve,
    template,
    ratio,
    leg,
    latex,
    unity,
  ]


def plot_comparison(
  run1_data: ROOT.TH1,
  run1_fit: dict,
  run2_data: ROOT.TH1,
  run2_fit: dict,
  rho: float,
  n_output_bins: int,
  outfile: str,
) -> None:
  canvas = ROOT.TCanvas(
    "c_kde_vs_mz_nominal",
    "KDE vs mz nominal",
    1500,
    650,
  )
  canvas.Divide(2, 1)

  keepalive: list = []
  keepalive.extend(
    _draw_run_panel(
      1,
      canvas,
      run1_data,
      run1_fit,
      f"Run 1 ({n_output_bins} bins)",
      rho,
    )
  )
  keepalive.extend(
    _draw_run_panel(
      2,
      canvas,
      run2_data,
      run2_fit,
      f"Run 2 ({n_output_bins} bins)",
      rho,
    )
  )
  canvas._keepalive = keepalive

  canvas.cd(0)
  canvas.Modified()
  canvas.Update()
  canvas.SaveAs(outfile)
  print(f"\nSaved {outfile}")


def output_png_path(n_output_bins: int) -> str:
  return os.path.join(
    OUTPUT_DIR,
    f"weighted_fixed_kde_vs_mz_nominal_{n_output_bins}bin.png",
  )


def main() -> int:
  if SOURCE_BINS % N_OUTPUT_BINS != 0:
    raise ValueError(
      f"N_OUTPUT_BINS={N_OUTPUT_BINS} must evenly divide SOURCE_BINS={SOURCE_BINS}"
    )

  rebin_factor = SOURCE_BINS // N_OUTPUT_BINS
  print(f"Combining {SOURCE_BINS} bins -> {N_OUTPUT_BINS} (factor {rebin_factor})")
  print(f"Loading KDE from {KDE_ROOT_FILE}")
  kde_shape, meta = load_kde_fit(KDE_ROOT_FILE)
  rho = meta["rho"]
  mix = meta.get("mix", float("nan"))
  train_alpha = meta.get("alpha", float("nan"))
  train_chi2 = meta.get("chi2", float("nan"))
  train_ndf = meta.get("ndf", float("nan"))
  print(
    f"Training fit: rho={rho:.6g}, mix={mix:.6g}, alpha={train_alpha:.6g}, "
    f"chi2/ndf={train_chi2 / max(train_ndf, 1):.6g} (ndf={train_ndf:.0f})"
  )

  datasets = [
    ("run1", RUN1_FILE),
    ("run2", RUN2_FILE),
  ]

  rebinned: dict[str, ROOT.TH1] = {}
  all_sumw2_ok = True

  for label, path in datasets:
    print(f"\nReading {HIST_NAME!r} from {path}")
    raw = open_projected_histogram(path, HIST_NAME, PROJECTION_AXIS)
    print(f"  projected {PROJECTION_AXIS}: {raw.GetNbinsX()} bins")
    all_sumw2_ok &= verify_sumw2_errors(raw, f"{label} ({SOURCE_BINS} bins)")

    rebinned[label] = rebin_histogram(raw, N_OUTPUT_BINS, f"{label}_{N_OUTPUT_BINS}bin")
    all_sumw2_ok &= verify_sumw2_errors(
      rebinned[label], f"{label} ({N_OUTPUT_BINS} bins, after rebin)"
    )

  if not all_sumw2_ok:
    print("\nWarning: at least one Sumw2 error check failed.")

  run1_fit = evaluate_fit("run1", rebinned["run1"], kde_shape)
  run2_fit = evaluate_fit("run2", rebinned["run2"], kde_shape)

  outfile = output_png_path(N_OUTPUT_BINS)
  os.makedirs(OUTPUT_DIR, exist_ok=True)
  plot_comparison(
    rebinned["run1"],
    run1_fit,
    rebinned["run2"],
    run2_fit,
    rho,
    N_OUTPUT_BINS,
    outfile,
  )
  return 0


if __name__ == "__main__":
  sys.exit(main())
