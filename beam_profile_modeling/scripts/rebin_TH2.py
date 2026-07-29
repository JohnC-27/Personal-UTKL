#!/usr/bin/env python3
"""
Rebin the high-resolution TH2s histograms and
write lower-resolution copies.
"""

import math
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Target number of bins per axis after combining the 2000-bin histograms.
N_OUTPUT_BINS = 100

# Assumes square TH2
SOURCE_BINS = 200
ROOT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "root_files")

SOURCE_FILES = {
  1: os.path.join(ROOT_FILES_DIR, "mz_nominal_2000bin_run1.root"),
  2: os.path.join(ROOT_FILES_DIR, "mz_nominal_2000bin_run2.root"),
}


def output_path(run: int, n_output_bins: int) -> str:
  return os.path.join(
    ROOT_FILES_DIR, f"mz_nominal_{n_output_bins}bin_run{run}.root"
  )


def list_th2_keys(filepath: str) -> list[str]:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  names: list[str] = []
  for key in tfile.GetListOfKeys():
    obj = key.ReadObj()
    if obj.InheritsFrom("TH2"):
      names.append(obj.GetName())
  tfile.Close()

  if not names:
    raise KeyError(f"no TH2 objects found in {filepath}")
  return sorted(names)


def open_histogram2d(filepath: str, hist_name: str) -> ROOT.TH2:
  tfile = ROOT.TFile.Open(filepath, "READ")
  if not tfile or tfile.IsZombie():
    raise OSError(f"cannot open ROOT file: {filepath}")

  obj = tfile.Get(hist_name)
  if not obj:
    tfile.Close()
    raise KeyError(f"object {hist_name!r} not found in {filepath}")

  if not obj.InheritsFrom("TH2"):
    tfile.Close()
    raise TypeError(f"{hist_name!r} is not TH2")

  hist = obj
  hist.SetDirectory(0)
  tfile.Close()
  return hist


def rebin_histogram2d(hist: ROOT.TH2, n_output_bins: int, name: str) -> ROOT.TH2:
  if SOURCE_BINS % n_output_bins != 0:
    raise ValueError(
      f"N_OUTPUT_BINS={n_output_bins} must evenly divide SOURCE_BINS={SOURCE_BINS}"
    )
  if hist.GetNbinsX() != SOURCE_BINS or hist.GetNbinsY() != SOURCE_BINS:
    raise ValueError(
      f"expected {SOURCE_BINS}x{SOURCE_BINS} bins before rebinning, got "
      f"{hist.GetNbinsX()}x{hist.GetNbinsY()}"
    )

  factor = SOURCE_BINS // n_output_bins
  rebinned = hist.Rebin2D(factor, factor, name)
  rebinned.SetDirectory(0)
  return rebinned


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


def rebin_run(run: int, n_output_bins: int) -> str:
  source_file = SOURCE_FILES[run]
  outfile = output_path(run, n_output_bins)
  hist_names = list_th2_keys(source_file)

  rebinned_hists: list[ROOT.TH2] = []
  for hist_name in hist_names:
    hist = open_histogram2d(source_file, hist_name)
    rebinned = rebin_histogram2d(hist, n_output_bins, hist_name)
    verify_sumw2_errors(rebinned, f"run{run} {hist_name}")
    rebinned_hists.append(rebinned)

  fout = ROOT.TFile.Open(outfile, "RECREATE")
  if not fout or fout.IsZombie():
    raise OSError(f"cannot create output file: {outfile}")

  for hist in rebinned_hists:
    hist.SetDirectory(fout)
    hist.Write()

  fout.Write()
  fout.Close()
  return outfile


def main() -> int:
  print(
    f"Rebinning {SOURCE_BINS}x{SOURCE_BINS} TH2s to "
    f"{N_OUTPUT_BINS}x{N_OUTPUT_BINS} bins"
  )

  for run in sorted(SOURCE_FILES):
    outfile = rebin_run(run, N_OUTPUT_BINS)
    print(f"Wrote {outfile}")

  return 0


if __name__ == "__main__":
  sys.exit(main())
