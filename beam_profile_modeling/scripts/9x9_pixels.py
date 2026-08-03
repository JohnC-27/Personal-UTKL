"""
Make an array of positions that the pixels are 
"""

import math
import numpy as np

pixel_spacing = 165.1   # mm
total_width = 101.6     # mm
detecting_width = 76.2  # mm
det_half_width = detecting_width / 2  # mm


positive_xpos = [pixel_spacing*i for i in range(0, 4)]
positive_ypos = [pixel_spacing*i for i in range(0, 4)]

print(positive_xpos)