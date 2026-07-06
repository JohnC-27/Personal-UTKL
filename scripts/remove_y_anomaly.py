"""
Resize TH2 histograms by trimming bins outside requested axis bounds.

By default, nominal beam-position histograms are trimmed at y >= 80 cm to remove
suspected non-physical beam behavior. Bounds for both x and y axes can be set
independently via resize_th2().
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

OUTPUT_ROOT_FILE1 = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "nominal_75x75.root"
)
OUTPUT_ROOT_FILE2 = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "mz_nominal_2000bin_run1_75x75.root"
)
OUTPUT_ROOT_FILE3 = os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "root_files", "mz_nominal_2000bin_run2_75x75.root"
)






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


def trim_hist_y_max(hist: ROOT.TH2, y_max_cm: float = Y_MAX_CM) -> ROOT.TH2:
  return resize_th2(hist, y_max=y_max_cm)


def write_resized_hists(
  hists: list[ROOT.TH2],
  output_path: str,
  x_min: float | None = None,
  x_max: float | None = None,
  y_min: float | None = None,
  y_max: float | None = None,
) -> None:
  outfile = ROOT.TFile.Open(output_path, "RECREATE")
  if not outfile or outfile.IsZombie():
    raise OSError(f"cannot create {output_path}")

  for hist in hists:
    resize_th2(hist, x_min=-75, x_max=75, y_min=-75, y_max=75).Write()

  outfile.Close()
  print(f"Wrote {len(hists)} histograms to {output_path}")


def write_trimmed_hists(hists: list[ROOT.TH2], output_path: str) -> None:
  write_resized_hists(hists, output_path, y_max=Y_MAX_CM)


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

  write_trimmed_hists(nominal_hists, OUTPUT_ROOT_FILE1)
  write_trimmed_hists(run1_hists, OUTPUT_ROOT_FILE2)
  write_trimmed_hists(run2_hists, OUTPUT_ROOT_FILE3)
  return 0


if __name__ == "__main__":
  sys.exit(main())
