import numpy as np
import scipy
import matplotlib.pyplot as plt
import os

def split_channels_numpy(csv_file):
  """
  Reads a CSV file with rows of 96 elements and splits each row
  into 4 sub-arrays of 24 elements.

  Returns a dictionary of NumPy arrays with shape (N_rows, 24).
  """
  skip = {48}
  data = np.genfromtxt(
  csv_file,
  delimiter=",",
  dtype=float,
  usecols=[i for i in range(97) if i not in skip],
  missing_values='',
  filling_values=np.nan,
)

  if data.shape[1] != 96:
    raise ValueError(f"Expected 96 columns, got {data.shape[1]}")

  channel_sets = {
    "set_1": data[:, 0:24],
    "set_2": data[:, 24:48],
    "set_3": data[:, 48:72],
    "set_4": data[:, 72:96],
  }

  return channel_sets

def get_mins(channels):
  mins = []
  for key in ['set_1', 'set_2', 'set_3', 'set_4']:
    data = channels[key]
    idx = np.argmax(np.abs(data), axis=1)
    row_indices = np.arange(data.shape[0])
    mins.append(data[row_indices, idx])
  return mins

def get_max(channels):
  maxs = [np.max(channels['set_1'],axis=1), np.max(channels['set_2'],axis=1), np.max(channels['set_3'],axis=1), np.max(channels['set_4'],axis=1)]
  return maxs