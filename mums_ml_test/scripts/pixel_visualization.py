"""
pixel_visualization.py — 2D x-y map of the pixel array with sense pads filled in.

Geometry (cm) comes from the pixel_params.csv row matching PIXEL_LABEL. Pixels are
centered on (0, 0) at the center-to-center spacing; each sense pad is a filled square.
If the row has a numeric pixel width, the full pixel outline is drawn around each pad.
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

# ===============================SCRIPT PARAMS===============================
PIXEL_LABEL = "4.0in"  # 2.5in, 3.0in, 4.0in
PIXEL_PARAMS_FILE = os.path.join(os.path.dirname(__file__), "..", "pixel_params.csv")
PLOT_DIR = os.path.join(os.path.dirname(__file__), "..", "plots")
# ===========================================================================


def read_pixel_params(filepath: str, label: str) -> dict:
  """Return spacing (cm), sense-pad width (cm), pixel width (cm), and pixels per side."""
  with open(filepath, newline="") as f:
    for row in csv.DictReader(f, skipinitialspace=True):
      if row["Label"] == label:
        return {
          "spacing": float(row["Center-to-center spacing cm"]),
          "sense_width": float(row["Sense-pad width cm"]),
          "pixel_width": float(row["Pixel width cm"]),
          "n_pixels": int(row["Num pixels"]),
        }
  raise ValueError(f"no row with Label {label!r} in {filepath}")


def plot_pixel_map(params: dict, label: str, out_path: str) -> None:
  spacing = params["spacing"]
  sense = params["sense_width"]
  pixel = params["pixel_width"]
  n_pixels = params["n_pixels"]
  centers = (np.arange(n_pixels) - (n_pixels - 1) / 2.0) * spacing

  fig, ax = plt.subplots(figsize=(7, 7))
  for cy in centers:
    for cx in centers:
      if pixel is not None:
        ax.add_patch(Rectangle(
          (cx - pixel / 2.0, cy - pixel / 2.0), pixel, pixel,
          fill=False, edgecolor="black", linewidth=0.8,
        ))
      ax.add_patch(Rectangle(
        (cx - sense / 2.0, cy - sense / 2.0), sense, sense,
        facecolor="tab:orange", edgecolor="black", linewidth=0.5,
      ))

  half = centers[-1] + max(sense, pixel or 0.0) / 2.0 + spacing / 4.0
  ax.set_xlim(-half, half)
  ax.set_ylim(-half, half)
  ax.set_aspect("equal")
  ax.set_xticks(centers)
  ax.set_yticks(centers)
  ax.tick_params(axis="x", labelrotation=45)
  ax.set_xlabel("x [cm]")
  ax.set_ylabel("y [cm]")
  pixel_text = f"{pixel:g} cm" if pixel is not None else "n/a"
  ax.set_title(
    f"{label}: {n_pixels}x{n_pixels}, spacing {spacing:g} cm, "
    f"sense pad {sense:g} cm, pixel {pixel_text}"
  )
  fig.tight_layout()
  fig.savefig(out_path, dpi=150)
  plt.close(fig)


def main() -> None:
  params = read_pixel_params(PIXEL_PARAMS_FILE, PIXEL_LABEL)
  os.makedirs(PLOT_DIR, exist_ok=True)  
  out_path = os.path.join(PLOT_DIR, f"pixel_map_{PIXEL_LABEL.replace('.','p')}.pdf")
  plot_pixel_map(params, PIXEL_LABEL, out_path)
  print(f"wrote {out_path}")


if __name__ == "__main__":
  main()
