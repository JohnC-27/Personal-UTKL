#!/usr/bin/env python3
"""
Shift and/or linearly combine KDE TH2s from 2d_kde fit ROOT files.

Modes (MODE):
  "shift"   — translate kde_shape / kde_template by (SHIFT_X, SHIFT_Y).
  "combine" — shift KDE_A and KDE_B independently, then form
              c1*shift(A) + c2*shift(B) as a standalone fit product.

Both modes write a full 2d_kde-layout ROOT (same keys) so plot scripts can
load the result as a normal fit product. Scan graphs and base metadata are
copied from the primary input; only the KDE TH2s are rebuilt.
"""

import os
import sys
from dataclasses import dataclass

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


# ===============================SCRIPT PARAMS===============================
# "shift" | "combine"
MODE = "combine"

INPUT_FILE_NAME = "jan2026_mm1_2d_kde"
INPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__), "..", "root_files", f"{INPUT_FILE_NAME}.root"
)

# Enter as mm. *0.1 converts to cm (nominal axis units).
SHIFT_X = 0.001 * 0.1
SHIFT_Y = 0.0 * 0.1

# Combine mode: load two fit products, shift each, then c1*A' + c2*B'.
# Primary (A) supplies target_hist, scan graphs, and base fit_meta.
COMBINE_FILE_A_NAME = "ml_nominal_mm1_2d_kde"  #.root added below
COMBINE_FILE_B_NAME = "pX_100um_mm1_2d_kde"
COMBINE_ROOT_A = os.path.join(
  os.path.dirname(__file__), "..", "root_files", f"{COMBINE_FILE_A_NAME}.root"
)
COMBINE_ROOT_B = os.path.join(
  os.path.dirname(__file__), "..", "root_files", f"{COMBINE_FILE_B_NAME}.root"
)

SHIFT_A_X = -0.0912 
SHIFT_A_Y = 0.0 
SHIFT_B_X = 0.06108 
SHIFT_B_Y = 0.0
COEFF_A = 0.5
COEFF_B = 0.5

OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__),
  "..",
  "root_files",
  f"{INPUT_FILE_NAME}_shifted_({SHIFT_X:g},{SHIFT_Y:g}).root",
)
COMBINE_OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(__file__),
  "..",
  "root_files",
  f"{COMBINE_FILE_A_NAME}_plus_{COMBINE_FILE_B_NAME}"
  f"_sA({SHIFT_A_X:g},{SHIFT_A_Y:g})"
  f"_sB({SHIFT_B_X:g},{SHIFT_B_Y:g})"
  f"_c({COEFF_A:g},{COEFF_B:g}).root",
)

# TH2s to translate / combine; written under the same names in the output file.
SHIFT_HIST_NAMES = ("kde_shape", "kde_template")
# ===========================================================================

REQUIRED_KEYS = (
  "target_hist",
  "chi2_vs_rho",
  "reduced_chi2_vs_rho",
  "alpha_vs_rho",
  "kde_shape",
  "kde_template",
  "fit_meta",
  "stats_meta",
)


@dataclass
class Th2Stats:
  integral: float
  mean_x: float
  mean_y: float


def in_axis_range(axis: ROOT.TAxis, value: float) -> bool:
  return axis.GetXmin() <= value <= axis.GetXmax()


def shift_th2(
  hist: ROOT.TH2,
  shift_x: float,
  shift_y: float,
  out_name: str,
) -> ROOT.TH2:
  """Translate content by (+shift_x, +shift_y); same binning as input."""
  out = hist.Clone(out_name)
  out.SetDirectory(0)
  out.SetTitle(
    f"{hist.GetTitle()} shifted ({shift_x:g},{shift_y:g})"
  )

  xaxis = hist.GetXaxis()
  yaxis = hist.GetYaxis()
  for ix in range(1, hist.GetNbinsX() + 1):
    x = xaxis.GetBinCenter(ix)
    src_x = x - shift_x
    for iy in range(1, hist.GetNbinsY() + 1):
      y = yaxis.GetBinCenter(iy)
      src_y = y - shift_y
      if in_axis_range(xaxis, src_x) and in_axis_range(yaxis, src_y):
        out.SetBinContent(ix, iy, hist.Interpolate(src_x, src_y))
      else:
        out.SetBinContent(ix, iy, 0.0)

  out.ResetStats()
  return out


def require_matching_binning(a: ROOT.TH2, b: ROOT.TH2) -> None:
  if a.GetNbinsX() != b.GetNbinsX() or a.GetNbinsY() != b.GetNbinsY():
    raise ValueError(
      f"bin count mismatch: {a.GetName()} "
      f"({a.GetNbinsX()}x{a.GetNbinsY()}) vs {b.GetName()} "
      f"({b.GetNbinsX()}x{b.GetNbinsY()})"
    )
  ax, bx = a.GetXaxis(), b.GetXaxis()
  ay, by = a.GetYaxis(), b.GetYaxis()
  edges = (
    (ax.GetXmin(), bx.GetXmin(), "x min"),
    (ax.GetXmax(), bx.GetXmax(), "x max"),
    (ay.GetXmin(), by.GetXmin(), "y min"),
    (ay.GetXmax(), by.GetXmax(), "y max"),
  )
  for va, vb, label in edges:
    if abs(va - vb) > 1e-9 * max(1.0, abs(va), abs(vb)):
      raise ValueError(f"axis {label} mismatch: {va} vs {vb}")


def add_th2s(
  hist_a: ROOT.TH2,
  hist_b: ROOT.TH2,
  coeff_a: float,
  coeff_b: float,
  out_name: str,
) -> ROOT.TH2:
  """Return coeff_a * hist_a + coeff_b * hist_b on hist_a's grid."""
  require_matching_binning(hist_a, hist_b)
  out = hist_a.Clone(out_name)
  out.SetDirectory(0)
  out.SetTitle(
    f"{coeff_a:g}*{hist_a.GetName()} + {coeff_b:g}*{hist_b.GetName()}"
  )
  for ix in range(1, out.GetNbinsX() + 1):
    for iy in range(1, out.GetNbinsY() + 1):
      value = (
        coeff_a * hist_a.GetBinContent(ix, iy)
        + coeff_b * hist_b.GetBinContent(ix, iy)
      )
      out.SetBinContent(ix, iy, value)
  out.ResetStats()
  return out


def th2_distribution_stats(hist: ROOT.TH2) -> Th2Stats:
  integral = 0.0
  mean_x_num = 0.0
  mean_y_num = 0.0
  for ix in range(1, hist.GetNbinsX() + 1):
    x = hist.GetXaxis().GetBinCenter(ix)
    for iy in range(1, hist.GetNbinsY() + 1):
      weight = hist.GetBinContent(ix, iy)
      if weight <= 0:
        continue
      y = hist.GetYaxis().GetBinCenter(iy)
      integral += weight
      mean_x_num += weight * x
      mean_y_num += weight * y

  if integral <= 0:
    return Th2Stats(0.0, 0.0, 0.0)
  return Th2Stats(
    integral=integral,
    mean_x=mean_x_num / integral,
    mean_y=mean_y_num / integral,
  )


def stats_meta_string(hist_stats: Th2Stats, kde_stats: Th2Stats) -> str:
  return (
    f"hist_integral={hist_stats.integral};"
    f"hist_mean_x={hist_stats.mean_x};hist_mean_y={hist_stats.mean_y};"
    f"kde_integral={kde_stats.integral};"
    f"kde_mean_x={kde_stats.mean_x};kde_mean_y={kde_stats.mean_y}"
  )


def load_fit_objects(filepath: str) -> dict[str, ROOT.TObject]:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  objects: dict[str, ROOT.TObject] = {}
  for name in REQUIRED_KEYS:
    obj = tfile.Get(name)
    if not obj:
      tfile.Close()
      raise KeyError(f"missing {name!r} in {filepath}")
    clone = obj.Clone(name)
    if hasattr(clone, "SetDirectory"):
      clone.SetDirectory(0)
    objects[name] = clone

  tfile.Close()
  return objects


def shift_named_hists(
  objects: dict[str, ROOT.TObject],
  shift_x: float,
  shift_y: float,
) -> dict[str, ROOT.TH2]:
  shifted: dict[str, ROOT.TH2] = {}
  for name in SHIFT_HIST_NAMES:
    hist = objects[name]
    if not hist.InheritsFrom("TH2"):
      raise TypeError(f"{name!r} is not TH2")
    shifted[name] = shift_th2(hist, shift_x, shift_y, name)
  return shifted


def augment_fit_meta_shift(
  meta: ROOT.TNamed,
  shift_x: float,
  shift_y: float,
) -> ROOT.TNamed:
  title = meta.GetTitle()
  for key in ("shift_x=", "shift_y="):
    if key in title:
      raise ValueError(f"fit_meta already contains {key!r}; refusing to double-shift")
  return ROOT.TNamed(
    "fit_meta",
    f"{title};shift_x={shift_x};shift_y={shift_y}",
  )


def augment_fit_meta_combine(
  meta: ROOT.TNamed,
  file_a: str,
  file_b: str,
  shift_a_x: float,
  shift_a_y: float,
  shift_b_x: float,
  shift_b_y: float,
  coeff_a: float,
  coeff_b: float,
) -> ROOT.TNamed:
  title = meta.GetTitle()
  for key in ("combine_coeff_a=", "combine_coeff_b="):
    if key in title:
      raise ValueError(f"fit_meta already contains {key!r}; refusing to re-combine")
  base_a = os.path.basename(file_a)
  base_b = os.path.basename(file_b)
  return ROOT.TNamed(
    "fit_meta",
    f"{title};combine_file_a={base_a};combine_file_b={base_b};"
    f"shift_a_x={shift_a_x};shift_a_y={shift_a_y};"
    f"shift_b_x={shift_b_x};shift_b_y={shift_b_y};"
    f"combine_coeff_a={coeff_a};combine_coeff_b={coeff_b}",
  )


def write_fit_product(
  outfile: str,
  objects: dict[str, ROOT.TObject],
  kde_hists: dict[str, ROOT.TH2],
  fit_meta: ROOT.TNamed,
  stats_meta: ROOT.TNamed,
) -> None:
  out_dir = os.path.dirname(os.path.abspath(outfile))
  os.makedirs(out_dir, exist_ok=True)

  fout = ROOT.TFile.Open(outfile, "RECREATE")
  if not fout or fout.IsZombie():
    raise OSError(f"cannot create output file: {outfile}")

  fout.cd()
  for name in ("target_hist", "chi2_vs_rho", "reduced_chi2_vs_rho", "alpha_vs_rho"):
    obj = objects[name].Clone(name)
    if hasattr(obj, "SetDirectory"):
      obj.SetDirectory(fout)
    obj.Write()

  for name, hist in kde_hists.items():
    hist.SetName(name)
    hist.SetDirectory(fout)
    hist.Write()

  fit_meta.Write()
  stats_meta.Write()
  fout.Write()
  fout.Close()


def run_shift() -> int:
  print(
    f"Mode: shift\n"
    f"Reading fit product from {INPUT_ROOT_FILE}\n"
    f"Shift (x, y) = ({SHIFT_X:g}, {SHIFT_Y:g}) cm\n"
    f"Writing {OUTPUT_ROOT_FILE}"
  )

  objects = load_fit_objects(INPUT_ROOT_FILE)
  shifted = shift_named_hists(objects, SHIFT_X, SHIFT_Y)

  target = objects["target_hist"]
  hist_stats = th2_distribution_stats(target)
  kde_stats = th2_distribution_stats(shifted["kde_template"])
  fit_meta = augment_fit_meta_shift(objects["fit_meta"], SHIFT_X, SHIFT_Y)
  stats_meta = ROOT.TNamed("stats_meta", stats_meta_string(hist_stats, kde_stats))

  write_fit_product(OUTPUT_ROOT_FILE, objects, shifted, fit_meta, stats_meta)

  print(
    f"Wrote {OUTPUT_ROOT_FILE}\n"
    f"  kde_template integral={kde_stats.integral:.6g} "
    f"mean=({kde_stats.mean_x:.6g}, {kde_stats.mean_y:.6g})"
  )
  return 0


def run_combine() -> int:
  print(
    f"Mode: combine (shift then sum)\n"
    f"A: {COMBINE_ROOT_A}\n"
    f"   shift=({SHIFT_A_X:g}, {SHIFT_A_Y:g}) cm  coeff={COEFF_A:g}\n"
    f"B: {COMBINE_ROOT_B}\n"
    f"   shift=({SHIFT_B_X:g}, {SHIFT_B_Y:g}) cm  coeff={COEFF_B:g}\n"
    f"Writing {COMBINE_OUTPUT_ROOT_FILE}"
  )

  objects_a = load_fit_objects(COMBINE_ROOT_A)
  objects_b = load_fit_objects(COMBINE_ROOT_B)
  shifted_a = shift_named_hists(objects_a, SHIFT_A_X, SHIFT_A_Y)
  shifted_b = shift_named_hists(objects_b, SHIFT_B_X, SHIFT_B_Y)

  combined: dict[str, ROOT.TH2] = {}
  for name in SHIFT_HIST_NAMES:
    combined[name] = add_th2s(
      shifted_a[name], shifted_b[name], COEFF_A, COEFF_B, name
    )

  target = objects_a["target_hist"]
  hist_stats = th2_distribution_stats(target)
  kde_stats = th2_distribution_stats(combined["kde_template"])
  fit_meta = augment_fit_meta_combine(
    objects_a["fit_meta"],
    COMBINE_ROOT_A,
    COMBINE_ROOT_B,
    SHIFT_A_X,
    SHIFT_A_Y,
    SHIFT_B_X,
    SHIFT_B_Y,
    COEFF_A,
    COEFF_B,
  )
  stats_meta = ROOT.TNamed("stats_meta", stats_meta_string(hist_stats, kde_stats))

  write_fit_product(
    COMBINE_OUTPUT_ROOT_FILE, objects_a, combined, fit_meta, stats_meta
  )

  print(
    f"Wrote {COMBINE_OUTPUT_ROOT_FILE}\n"
    f"  kde_template integral={kde_stats.integral:.6g} "
    f"mean=({kde_stats.mean_x:.6g}, {kde_stats.mean_y:.6g})"
  )
  return 0


def main() -> int:
  if MODE == "shift":
    return run_shift()
  if MODE == "combine":
    return run_combine()
  raise ValueError(f"MODE must be 'shift' or 'combine', got {MODE!r}")


if __name__ == "__main__":
  sys.exit(main())
