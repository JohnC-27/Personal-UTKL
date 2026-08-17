import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import sys

sys.path.append('..')

from readout_functions import *

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")


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
  old_data_dir = os.path.join(os.path.dirname(__file__), "..", "2026-07-07")
  data_dir = os.path.join(os.path.dirname(__file__), "..", "2026-07-08")

  filelist1 = []
  for files in os.listdir(old_data_dir):
    if files.endswith(".csv"):
      filelist1.append(os.path.join(old_data_dir, files))

  filelist = []
  for files in os.listdir(data_dir):
    if files.endswith(".csv"):
      filelist.append(os.path.join(data_dir, files))

  filelist.sort(key = lambda x: int(x.split('_')[1].split('.')[0]))  # Sort by the number in the filename

  filelist1.sort(key = lambda x: int(x.split('_')[1].split('.')[0]))

  int_times = np.array([1,2,3,4,5,6,7,8])
  int_times_2 = np.array([10, 20, 40, 60, 80, 100, 120, 140])
  int_times_3 = np.array([26, 28, 30, 32, 34, 36, 38, 40])


  fig, axs = plt.subplots(figsize=(10,8))

  #print(split_channels_numpy(filelist[0]))
  mins = [get_mins(split_channels_numpy(file))[1] for file in filelist] #[np.min(channels['set_1'],axis=1), np.min(channels['set_2'],axis=1), np.min(channels['set_3'],axis=1), np.min(channels['set_4'],axis=1)]
  maxes = [get_max(split_channels_numpy(file))[1] for file in filelist]

  mins2 = [get_mins(split_channels_numpy(file))[1] for file in filelist1]

  mins_1000 = [get_mins(split_channels_numpy(file))[3] for file in filelist]

  popt, pcov = fit_to_data(int_times_2, mins[1])
  print(popt[0]/1.e-6*100.e-12)

  popt1, pcov1 = fit_to_data(int_times_2, mins[2])
  print(popt1[0]/1.e-6*100.e-12)

  popt2, pcov2 = fit_to_data(int_times_2[1:], mins[5][1:])
  print(popt2[0]/1.e-6*100.e-12)

  popt3, pcov3 = fit_to_data(int_times_2[1:], mins[5][1:])
  print(popt3[0]/1.e-6*100.e-12)

  #stds = [np.std(channels['set_1'],axis=1), np.std(channels['set_2'],axis=1), np.std(channels['set_3'],axis=1), np.std(channels['set_4'],axis=1)]

    

  '''axs.plot(int_times_2, mins2[4], marker='o', label="100 us, old setup")
  axs.plot(int_times_2, mins[12], marker='o', label="100 us, new setup")
  axs.plot(int_times_2, mins2[5], marker='o', label="20 us, old setup")
  axs.plot(int_times_2, mins[1], marker='o', label="20 us, new setup")

  axs.grid()

  axs.legend(fontsize=14)
  axs.set_xlabel("Integration time (us)")
  axs.set_ylabel("SWIC scanner signal (V)")
  fig.suptitle("Integration time scans, old vs new setups")

  fig.savefig("Figure_1.pdf")'''

  fig1, axs1 = plt.subplots(figsize=(10,8))

  axs1.plot(int_times*10., mins2[2], marker='o', label="1.4 ft RG-58")
  axs1.plot(int_times*10., mins2[3], marker='o', label="2.7 ft RG-58")
  axs1.plot(int_times*10., mins2[4], marker='o', label="5.3 ft RG-58")
  axs1.plot(int_times*10., mins2[5], marker='o', label="11 ft RG-58")
  axs1.plot(int_times*10., mins2[6], marker='o', label="21 ft RG-58")
  axs1.plot(int_times*10., mins[0], marker='o', label="43 ft RG-58")
  axs1.plot(int_times*10., mins[1], marker='o', label="67 ft RG-58")
  axs1.plot(int_times*10., mins[2], marker='o', label="160 ft RG-58")

  axs1.grid()

  axs1.legend(fontsize=14)
  axs1.set_xlabel("Integration time (us)")
  axs1.set_ylabel("SWIC scanner signal (V)")
  #fig1.suptitle("9 uA 30 us wide pulse, different cable lengths")

  fig1.savefig(os.path.join(OUTPUT_DIR, "cable_lengths_07-08.pdf"))

  fig2, axs2 = plt.subplots(figsize=(10,8))

  axs2.plot(int_times_3, mins2[8], marker='o', label="1.4 ft RG-58")
  #axs2.plot(int_times_2, linear_fit(int_times_2, popt[0], popt[1]), linestyle='--', label="1.15 uA best fit")
  axs2.plot(int_times_3, mins2[9], marker='o', label="2.7 ft RG-58")
  #axs2.plot(int_times_2, linear_fit(int_times_2, popt1[0], popt1[1]), linestyle='--', label="2.02 uA best fit")
  axs2.plot(int_times_3, mins2[10], marker='o', label="5.3 ft RG-58")
  axs2.plot(int_times_3, mins2[11], marker='o', label="11 ft RG-58")
  axs2.plot(int_times_3, mins2[13], marker='o', label="21 ft RG-58")
  axs2.plot(int_times_3, mins[6], marker='o', label="21 ft TP", linestyle="--", color='yellowgreen')
  axs2.plot(int_times_3, mins[5], marker='o', label="43 ft RG-58")
  axs2.plot(int_times_3, mins[7], marker='o', label="50 ft TP", linestyle="--", color='teal')
  axs2.plot(int_times_3, mins[4], marker='o', label="67 ft RG-58")
  axs2.plot(int_times_3, mins[3], marker='o', label="160 ft RG-58")
  

  #axs2.plot(int_times_2, mins2[1], marker='o', label="2 V, 40 mV low")

  axs2.grid()

  axs2.legend(fontsize=14)
  axs2.set_xlabel("Integration time (us)")
  axs2.set_ylabel("SWIC scanner signal (V)")
  #fig2.suptitle("9 uA 30 us wide pulses, different cable lengths")

  fig2.savefig(os.path.join(OUTPUT_DIR, "both_cables_lengths_26-40us.pdf"))
