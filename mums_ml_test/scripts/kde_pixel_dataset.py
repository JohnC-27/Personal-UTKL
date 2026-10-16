"""
kde_pixel_dataset.py — Integrate three KDEs (MM1, MM2, MM3) over a 9x9 pixel array.

Import and call pixel_dataset_row(); it returns one float32 row:
  [MM1 pixels (n^2), MM2 pixels (n^2), MM3 pixels (n^2), x_shift_cm, y_shift_cm]
All lengths are in cm. Spacing, sense-pad width, and n (pixels per side) come from
the pixel_params.csv row matching the chosen Label; integration covers the sense pad.
Pixel values are events / 1e9 POT. Pixel order within a detector is
index = iy * n + ix, with ix and iy increasing along +x and +y.
Integration and sums are fp64; only the returned row is fp32.
"""

import csv
import os
from dataclasses import dataclass

import numpy as np
import ROOT

# ===============================MODULE PARAMS===============================
PIXEL_PARAMS_FILE = os.path.join(os.path.dirname(__file__), "..", "pixel_params.csv")
DEFAULT_PIXEL_LABEL = "3.0in" #2.5in, 3.0in, 4.0in
N_SUB = 10  # midpoint sub-points per axis within each pixel
DEFAULT_POT = 1.0e9
# ===========================================================================


@dataclass
class PixelKde:
  """One detector's KDE as built by 2d_kde.py: alpha * (mix * u + (1 - mix) * m)."""
  pdf_unmirrored: ROOT.RooNDKeysPdf
  pdf_mirrored: ROOT.RooNDKeysPdf
  mix: float
  alpha: float
  ctx: object  # 2d_kde.KdeFitContext the pdfs were built from


def read_pixel_params(filepath: str, label: str) -> dict[str, float]:
  """Return spacing (cm), sense-pad width (cm), and pixels per side for one Label row."""
  with open(filepath, newline="") as f:
    for row in csv.DictReader(f, skipinitialspace=True):
      if row["Label"] == label:
        return {
          "spacing": float(row["Center-to-center spacing cm"]),
          "sense_width": float(row["Sense-pad width cm"]),
          "n_pixels": int(row["Num pixels"]),
        }
  raise ValueError(f"no row with Label {label!r} in {filepath}")


def kde_pixel_values(kde: PixelKde, params: dict[str, float]) -> np.ndarray:
  """Expected events per pixel (unscaled by POT), fp64, flattened iy-major."""
  target = kde.ctx.target
  x_var, y_var, argset = kde.ctx.x_var, kde.ctx.y_var, kde.ctx.argset
  spacing = params["spacing"]
  width = params["sense_width"]
  n_pixels = params["n_pixels"]

  # alpha * pdf is counts per histogram bin; divide by bin area for counts/cm^2.
  bin_area = target.GetXaxis().GetBinWidth(1) * target.GetYaxis().GetBinWidth(1)
  centers = (np.arange(n_pixels) - (n_pixels - 1) / 2.0) * spacing
  offsets = ((np.arange(N_SUB) + 0.5) / N_SUB - 0.5) * width
  sub_area = (width / N_SUB) ** 2

  # RooRealVar.setVal clamps out-of-range values, which would silently corrupt edge pixels.
  edge = centers[-1] + width / 2.0
  for axis in (target.GetXaxis(), target.GetYaxis()):
    if axis.GetXmin() > -edge or axis.GetXmax() < edge:
      raise ValueError(f"pixel array (+/-{edge} cm) exceeds histogram range")

  pixels = np.zeros((n_pixels, n_pixels), dtype=np.float64)
  for iy, cy in enumerate(centers):
    for ix, cx in enumerate(centers):
      total = 0.0
      for y in cy + offsets:
        y_var.setVal(y)
        for x in cx + offsets:
          x_var.setVal(x)
          u = kde.pdf_unmirrored.getVal(argset)
          m = kde.pdf_mirrored.getVal(argset)
          total += kde.mix * u + (1.0 - kde.mix) * m
      pixels[iy, ix] = kde.alpha * total * sub_area / bin_area
  return pixels.flatten()


def pixel_dataset_row(
  kde_mm1: PixelKde,
  kde_mm2: PixelKde,
  kde_mm3: PixelKde,
  x_shift_cm: float,
  y_shift_cm: float,
  pot: float = DEFAULT_POT,
  pixel_label: str = DEFAULT_PIXEL_LABEL,
  params_file: str = PIXEL_PARAMS_FILE,
) -> np.ndarray:
  """Pixel integrals for MM1-MM3 in events / 1e9 POT, then the shift labels (cm), as fp32."""
  params = read_pixel_params(params_file, pixel_label)
  scale = 1.0e9 / pot
  parts = [kde_pixel_values(kde, params) * scale for kde in (kde_mm1, kde_mm2, kde_mm3)]
  parts.append(np.array([x_shift_cm, y_shift_cm], dtype=np.float64))
  return np.concatenate(parts).astype(np.float32)
