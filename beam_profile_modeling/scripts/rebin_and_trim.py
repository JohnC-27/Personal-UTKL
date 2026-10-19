#!/usr/bin/env python3
"""
Combines rebin_TH2.py and remove_y_anomaly.py.
Takes in "raw" TH2s from simulation
Rebin square TH2s to N_OUTPUT_BINS x N_OUTPUT_BINS, then trim the y axis
symmetrically to [-Y_MAX_CM, Y_MAX_CM] to remove the tracking plane
anomaly at the +y edge.
 
File-wise over a list of files. For each file, name the TH2s to process 
directly or set MUMS_ONLY = True to only process MuMs TH2s. Output is root files
in the same directory with the rebinned and resized TH2s.
Output naming: <stem>[_mums]_[binx]x[biny]_corrected.root in ROOT_FILES_DIR, where
<stem> is the matching OUTPUT_FILE_STEMS entry or the input file stem.
"""

import math
import os
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Target number of bins per axis after rebinning (must divide the source bins).
# Check source bins in the root file catalog or witht the inspect script
N_OUTPUT_BINS = 100
# Keep bins with centers in [-Y_MAX_CM, Y_MAX_CM]; x axis is not trimmed.
Y_MAX_CM = 95.0

# input and output directory
ROOT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "root_files")

# with .root
INPUT_FILE_NAMES = [
  "nominal_1e9pot_zarko.root",
  "protonbeam_pY_1sigma_1e9pot_zarko.root",
  "protonbeam_pY_2sigma_1e9pot_zarko.root",
  #"jan2026studies_nominal.root",
]
INPUT_FILES = [os.path.join(ROOT_FILES_DIR, f) for f in INPUT_FILE_NAMES]

# Output file name is <stem>[_mums]_{binx}x{biny}_corrected.root.
# An empty list or an empty entry uses the input file stem; otherwise the
# entry replaces it. Must keep in same order as INPUT_FILE_NAMES. 
OUTPUT_FILE_STEMS = [
  "nom_1e9pot_zarko",
  "pbeam_pY_1sigma_1e9pot_zarko",
  "pbeam_pY_2sigma_1e9pot_zarko",
] # list of strings or empty

# ---------------------------------------------------------------------------
# TH2 selection: key = input file basename, value = list of TH2 names.
# Files not listed here (or value = None) transform every TH2 in the file.
#
# Set MUMS_ONLY to only process MuMS hists. The flag overrides selections in
# HISTS_TO_TRANSFORM and None values
# ---------------------------------------------------------------------------
MUMS_ONLY = True
HISTS_TO_TRANSFORM: dict[str, list[str] | None] = {
  "nominal_1e9pot_zarko.root": [
    "Beamline_shift_pXxyposMM1",
    "Beamline_shift_pXxyposMM2",
    "Beamline_shift_pXxyposMM3",
  ],
  #"jan2026studies_nominal.root": ["nominal_xypos_1", "nominal_xypos_2", "nominal_xypos_3"],
}


def list_th2_keys(filepath: str, mums_flag: bool) -> list[str]:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  names: list[str] = []
  keys = tfile.GetListOfKeys()

  if mums_flag == True:
    for key in keys:
      obj = key.ReadObj()
      if obj.InheritsFrom("TH2") and 'alcove' in key.GetTitle().lower():
        names.append(key.GetName())
  else:
    for key in keys:
      obj = key.ReadObj()
      if obj.InheritsFrom("TH2"):
        names.append(key.GetName())
  tfile.Close()

  if not names:
    raise KeyError(f"no TH2 objects found in {filepath}")
  return sorted(names)


def selected_hist_names(filepath: str, mums_flag: bool) -> list[str]:
  selection = HISTS_TO_TRANSFORM.get(os.path.basename(filepath))
  if selection is None or mums_flag == True:
    return list_th2_keys(filepath, mums_flag)
  return list(selection)


def load_th2_group(filepath: str, hist_names: list[str]) -> list[ROOT.TH2]:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open {filepath}")

  hists: list[ROOT.TH2] = []
  for hist_name in hist_names:
    hist = tfile.Get(hist_name)
    if not hist or not hist.InheritsFrom("TH2"):
      tfile.Close()
      raise KeyError(f"missing or invalid TH2 {hist_name!r} in {filepath}")
    hist.SetDirectory(0)
    hists.append(hist)

  tfile.Close()
  return hists


def rebin_histogram2d(hist: ROOT.TH2, n_output_bins: int) -> ROOT.TH2:
  nx = hist.GetNbinsX()
  ny = hist.GetNbinsY()
  if nx != ny:
    raise ValueError(f"{hist.GetName()!r} is not square: {nx}x{ny} bins")
  if nx % n_output_bins != 0:
    raise ValueError(
      f"N_OUTPUT_BINS={n_output_bins} must evenly divide {nx} bins in {hist.GetName()!r}"
    )

  factor = nx // n_output_bins
  rebinned = hist.Rebin2D(factor, factor, hist.GetName())
  rebinned.SetDirectory(0)
  return rebinned


def _axis_bin_range(
  axis: ROOT.TAxis,
  low: float | None,
  high: float | None,
  axis_name: str,
  hist_name: str,
) -> tuple[int, int]:
  if low is None:
    low = axis.GetXmin()
  if high is None:
    high = axis.GetXmax()
  if low > high:
    raise ValueError(
      f"{axis_name} lower bound {low} exceeds upper bound {high} in {hist_name!r}"
    )

  min_i: int | None = None
  max_i = 0
  for i in range(1, axis.GetNbins() + 1):
    center = axis.GetBinCenter(i)
    if center >= low and min_i is None:
      min_i = i
    if center <= high:
      max_i = i

  if min_i is None or max_i == 0 or min_i > max_i:
    raise ValueError(
      f"no {axis_name} bins with centers in [{low}, {high}] in {hist_name!r}"
    )
  return min_i, max_i


def resize_th2(
  hist: ROOT.TH2,
  x_min: float | None = None,
  x_max: float | None = None,
  y_min: float | None = None,
  y_max: float | None = None,
) -> ROOT.TH2:
  """Return a copy of hist trimmed to the requested axis bounds."""
  hist_name = hist.GetName()
  xaxis = hist.GetXaxis()
  yaxis = hist.GetYaxis()

  min_ix, max_ix = _axis_bin_range(xaxis, x_min, x_max, "x", hist_name)
  min_iy, max_iy = _axis_bin_range(yaxis, y_min, y_max, "y", hist_name)

  out = hist.Clone(hist_name)
  out.SetDirectory(0)
  out.SetBins(
    max_ix - min_ix + 1,
    xaxis.GetBinLowEdge(min_ix),
    xaxis.GetBinUpEdge(max_ix),
    max_iy - min_iy + 1,
    yaxis.GetBinLowEdge(min_iy),
    yaxis.GetBinUpEdge(max_iy),
  )

  for ix in range(min_ix, max_ix + 1):
    for iy in range(min_iy, max_iy + 1):
      dest_ix = ix - min_ix + 1
      dest_iy = iy - min_iy + 1
      out.SetBinContent(dest_ix, dest_iy, hist.GetBinContent(ix, iy))
      out.SetBinError(dest_ix, dest_iy, hist.GetBinError(ix, iy))

  return out


def verify_sumw2_errors(hist: ROOT.TH2, label: str) -> bool:
  """Check that GetBinError(ix,iy) equals sqrt(sumw2[bin]) for every bin."""
  if hist.GetSumw2N() == 0:
    print(f"{label}: no Sumw2 array (GetSumw2N() == 0)")
    return False

  sw2 = hist.GetSumw2()
  max_abs_diff = 0.0
  n_mismatch = 0

  for ix in range(1, hist.GetNbinsX() + 1):
    for iy in range(1, hist.GetNbinsY() + 1):
      err = hist.GetBinError(ix, iy)
      expected = math.sqrt(max(sw2.At(hist.GetBin(ix, iy)), 0.0))
      abs_diff = abs(err - expected)
      if abs_diff > 1e-6 * max(abs(err), abs(expected), 1.0):
        n_mismatch += 1
      max_abs_diff = max(max_abs_diff, abs_diff)

  ok = n_mismatch == 0
  status = "PASS" if ok else "FAIL"
  print(
    f"{label}: Sumw2 error check {status} "
    f"({hist.GetNbinsX()}x{hist.GetNbinsY()} bins, mismatches={n_mismatch}, "
    f"max |diff|={max_abs_diff:.3e})"
  )
  return ok


def output_path(input_path: str, hists: list[ROOT.TH2]) -> str:
  shapes = {(h.GetNbinsX(), h.GetNbinsY()) for h in hists}
  if len(shapes) != 1:
    raise ValueError(f"transformed TH2s in {input_path} have mixed shapes: {shapes}")
  nx, ny = shapes.pop()

  input_name = os.path.basename(input_path)
  i = INPUT_FILE_NAMES.index(input_name)
  custom_stem = OUTPUT_FILE_STEMS[i] if OUTPUT_FILE_STEMS else ""
  stem = custom_stem or os.path.splitext(input_name)[0]
  if MUMS_ONLY:
    stem += "_mums"

  return os.path.join(ROOT_FILES_DIR, f"{stem}_{nx}x{ny}_corrected.root")


def transform_file(input_path: str, n_output_bins: int, y_max: float) -> str:
  hists = load_th2_group(input_path, selected_hist_names(input_path, MUMS_ONLY))

  corrected: list[ROOT.TH2] = []
  for hist in hists:
    rebinned = rebin_histogram2d(hist, n_output_bins)
    trimmed = resize_th2(rebinned, y_min=-y_max, y_max=y_max)
    verify_sumw2_errors(trimmed, f"{os.path.basename(input_path)} {trimmed.GetName()}")
    corrected.append(trimmed)

  outfile = output_path(input_path, corrected)
  fout = ROOT.TFile.Open(outfile, "RECREATE")
  if not fout or fout.IsZombie():
    raise OSError(f"cannot create output file: {outfile}")
  for hist in corrected:
    hist.SetDirectory(fout)
    hist.Write()
  fout.Close()
  print(f"Wrote {len(corrected)} histograms to {outfile}")
  return outfile


def main() -> int:
  if OUTPUT_FILE_STEMS and len(OUTPUT_FILE_STEMS) != len(INPUT_FILE_NAMES):
    raise ValueError(
      f"OUTPUT_FILE_STEMS has {len(OUTPUT_FILE_STEMS)} entries but "
      f"INPUT_FILE_NAMES has {len(INPUT_FILE_NAMES)}"
    )

  print(
    f"Rebinning TH2s to {N_OUTPUT_BINS}x{N_OUTPUT_BINS} bins, "
    f"then trimming y to [{-Y_MAX_CM}, {Y_MAX_CM}] cm"
  )
  for input_path in INPUT_FILES:
    transform_file(input_path, N_OUTPUT_BINS, Y_MAX_CM)
  return 0


if __name__ == "__main__":
  sys.exit(main())
