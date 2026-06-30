"""
Remove bins that contain suspected non-physical beam behavior. All bins with y coord >= 80cm
"""
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

Y_MAX_CM = 80.0

INPUT_ROOT_FILE1 = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "nominal.root"
)
INPUT_ROOT_FILE2 = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "mz_nominal_2000bin_run1.root"
)
INPUT_ROOT_FILE3 = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "mz_nominal_2000bin_run2.root"
)

OUTPUT_ROOT_FILE = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "nominal_corrected.root"
)


def corrected_output_path(input_path: str) -> str:
  root, ext = os.path.splitext(input_path)
  return f"{root}_corrected{ext}"


def load_th2_group(
  filepath: str,
  hist_names: list[str],
) -> list[ROOT.TH2]:
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


def trim_hist_y_max(hist: ROOT.TH2, y_max_cm: float = Y_MAX_CM) -> ROOT.TH2:
  max_iy = 0
  for iy in range(1, hist.GetNbinsY() + 1):
    if hist.GetYaxis().GetBinCenter(iy) <= y_max_cm:
      max_iy = iy

  if max_iy == 0:
    raise ValueError(f"no y bins with center <= {y_max_cm} cm in {hist.GetName()!r}")

  out = hist.Clone(hist.GetName())
  out.SetDirectory(0)
  out.SetBins(
    hist.GetNbinsX(),
    hist.GetXaxis().GetXmin(),
    hist.GetXaxis().GetXmax(),
    max_iy,
    hist.GetYaxis().GetXmin(),
    hist.GetYaxis().GetBinUpEdge(max_iy),
  )
  return out


def write_trimmed_hists(hists: list[ROOT.TH2], output_path: str) -> None:
  outfile = ROOT.TFile.Open(output_path, "RECREATE")
  if not outfile or outfile.IsZombie():
    raise OSError(f"cannot create {output_path}")

  for hist in hists:
    trim_hist_y_max(hist).Write()

  outfile.Close()
  print(f"Wrote {len(hists)} histograms to {output_path}")


def main() -> int:
  nominal_hists = load_th2_group(
    INPUT_ROOT_FILE1,
    ["nominalxyposMM1", "nominalxyposMM2", "nominalxyposMM3"],
  )
  run1_hists = load_th2_group(
    INPUT_ROOT_FILE2,
    ["NominalxyposMM1", "NominalxyposMM2", "NominalxyposMM3"],
  )
  run2_hists = load_th2_group(
    INPUT_ROOT_FILE3,
    ["NominalxyposMM1", "NominalxyposMM2", "NominalxyposMM3"],
  )

  write_trimmed_hists(nominal_hists, OUTPUT_ROOT_FILE)
  write_trimmed_hists(run1_hists, corrected_output_path(INPUT_ROOT_FILE2))
  write_trimmed_hists(run2_hists, corrected_output_path(INPUT_ROOT_FILE3))
  return 0


if __name__ == "__main__":
  sys.exit(main())