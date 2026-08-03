"""
Gets bin positions that correspond to pixel positions. Only takes bins that are 100%
within the bounds of the detecting size of the pixel as to not inflate statistics.
Assumes size of pixels is maxed out at 3in/76.2mm (detecting size). Pixel spacing variable
for optimization later.
"""
import numpy as np


#TODO Finish. Not sure of the return format. Ranges might be best
# so there can be a two sets of nested for loops that iterate over a list of index ranges. 
def get_bin_indicies_9x9_pixels(
  pixel_spacing: float, 
  total_width: float,
  detecting_width: float
  ) -> list:
  # Trying rounding down because the fractional portions dont even add up to a full mm over the entire grid,
  # and not attempting to interpolate between bins anyway.
  # keeping the pixels 76*76 is more important
  pixel_spacing = 165   # mm, actually 165.1mm
  total_width = 101.6     # mm
  detecting_width = 76  # mm, actually 76.2mm
  detecting_half_width = detecting_width / 2  # mm


  x_centers = [round(pixel_spacing*i,1) for i in range(-4, 5)]
  y_centers = [round(pixel_spacing*i,1) for i in range(-4, 5)]

  x_lower = [round(pixel_spacing*i - detecting_half_width,1) for i in range(-4, 5)]
  y_lower = [round(pixel_spacing*i - detecting_half_width,1) for i in range(-4, 5)]

  x_ranges = [(round(x-detecting_half_width,1), round(x+detecting_half_width,1)) for x in x_centers]
  y_ranges = [(round(y-detecting_half_width,1), round(y+detecting_half_width,1)) for y in y_centers]

  x_bins_pos = []
  y_bins_pos = []
  for x in x_lower:
    for i in range(0, 77):
      x_bins_pos.append(x+i)
      y_bins_pos.append(x+i)


  # converts from positinon in mm to bin index
  # FOR 2000 BINS over 2m x 2m.
  x_bins_indicies = [1000+i for i in x_bins_pos]
  y_bins_indicies = [1000+i for i in y_bins_pos]
  print(x_bins_indicies)
  return None


def get_cm_ranges_9x9_pixels(
  pixel_spacing: float = 16.51,
  detecting_width: float = 7.62,
) -> np.array:
  """INPUT IN CM. Returns ranges of pixels. Symmetric in x/y so returns one np.array. Rounded to .01cm"""

  detecting_half_width = detecting_width / 2  # cm

  x_centers = [round(pixel_spacing*i,2) for i in range(-4, 5)]
  x_ranges = [(round(x-detecting_half_width,2), round(x+detecting_half_width,2)) for x in x_centers]

  return np.array(x_ranges)


print(get_cm_ranges_9x9_pixels())