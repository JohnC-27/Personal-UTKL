import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import sys

sys.path.append('..')

from readout_functions import *

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
coax_data_dir = os.path.join(os.path.dirname(__file__), "..", "2026-07-13")
tp_data_dir = os.path.join(os.path.dirname(__file__), "..", "2026-08-11")


def linear_fit(x, m, b):
  return m*x + b

def fit_to_data(x_data, y_data):
  popt, pcov = scipy.optimize.curve_fit(linear_fit, x_data, y_data)
  return popt, pcov

plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 24

if __name__ == "__main__":

  linstyles_for_current = ['-','--',':']
  int_times = np.array([1,2,3,4,5,6,7,8])
  labels = ["1 mA", "2 mA", "5 mA", "10 mA"]

  coax_filelist = []
  for files in os.listdir(coax_data_dir):
    if files.endswith(".csv"):
      coax_filelist.append(os.path.join(coax_data_dir, files))

  tp_filelist = []
  for files in os.listdir(tp_data_dir):
    if files.endswith(".csv"):
      tp_filelist.append(os.path.join(tp_data_dir, files))


  coax_filelist.sort(key = lambda x: int(x.split('_')[1].split('.')[0]))
  tp_filelist.sort(key = lambda x: int(x.split('_')[1].split('.')[0]))

  # channel 1 is 100 pF board
  # 
  coax_mins = [get_mins(split_channels_numpy(file))[1] for file in coax_filelist]
  tp_mins = mins = [get_mins(split_channels_numpy(file))[1] for file in tp_filelist]

  

  fig4, axs4 = plt.subplots(figsize=(10,8))
  axs4.plot(int_times*10., coax_mins[4], marker='o', color=current_colors_rg_58[0],linestyle=linstyles_for_current[0],label="10 uA, 21 ft RG-58")
  axs4.plot(int_times*10., coax_mins[8], marker='o', color=current_colors_rg_58[0],linestyle=linstyles_for_current[1],label="10 uA, 65 ft RG-58")
  axs4.plot(int_times*10., tp_mins[0], marker='o', color=current_colors_tp[0],linestyle=linstyles_for_current[0],label="10 uA, 21 ft TP")
  axs4.plot(int_times*10., tp_mins[4], marker='o', color=current_colors_tp[0],linestyle=linstyles_for_current[1],label="10 uA, 50 ft TP")
  axs4.plot(int_times*10., coax_mins[5], marker='o', color=current_colors_rg_58[1],linestyle=linstyles_for_current[0],label="20 uA, 21 ft RG-58")
  axs4.plot(int_times*10., coax_mins[9], marker='o', color=current_colors_rg_58[1],linestyle=linstyles_for_current[1],label="20 uA, 65 ft RG-58")
  axs4.plot(int_times*10., tp_mins[1], marker='o', color=current_colors_tp[1],linestyle=linstyles_for_current[0],label="20 uA, 21 ft TP")
  axs4.plot(int_times*10., tp_mins[5], marker='o', color=current_colors_tp[1],linestyle=linstyles_for_current[1],label="20 uA, 50 ft TP")
  #axs4.plot(int_times*10., mins_100000[49], marker='o', color=colors_for_current[2],linestyle=linstyles_for_current[2])

  axs4.grid()

  axs4.legend(fontsize=14)
  axs4.set_xlabel("Integration time (us)")
  axs4.set_ylabel("SWIC scanner signal (V)")
  axs4.set_ylim(-6,0)
  #fig4.suptitle("30 us wide pulses into 100 pF board")

  fig4.savefig(os.path.join(OUTPUT_DIR, "100pF_rg-58_vs_tp.pdf"))