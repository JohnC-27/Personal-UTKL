import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import sys

sys.path.append('..')

from readout_functions import *


def linear_fit(x, m, b):
  return m*x + b

def fit_to_data(x_data, y_data):
  popt, pcov = scipy.optimize.curve_fit(linear_fit, x_data, y_data)
  return popt, pcov

plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 24


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
  filelist = []
  for file in os.listdir(SCRIPT_DIR):
    if file.endswith(".csv"):
      filelist.append(os.path.join(SCRIPT_DIR, file))

  filelist.sort(key=lambda x: int(os.path.basename(x).split('_')[1].split('.')[0]))


  # integration times
  int_times = np.array([10,20,30,40,50,60,70,80])
  #int_times_2 = np.array([10, 20, 40, 60, 80, 100, 120, 140])
  #int_times_3 = np.array([26, 28, 30, 32, 34, 36, 38, 40])


  fig, axs = plt.subplots(figsize=(10,8))
  print(filelist)
  #print(split_channels_numpy(filelist[0]))
  mins = [get_mins(split_channels_numpy(file))[1] for file in filelist] #[np.min(channels['set_1'],axis=1), np.min(channels['set_2'],axis=1), np.min(channels['set_3'],axis=1), np.min(channels['set_4'],axis=1)]
  maxes = [get_max(split_channels_numpy(file))[1] for file in filelist]

  mins_1000 = [get_mins(split_channels_numpy(file))[3] for file in filelist]

  # to cut out anomaly on the higher current scans
  fit_mask_10_50 = (int_times >= 10) & (int_times <= 50)
  fit_times_10_50 = int_times[fit_mask_10_50]

  popt, pcov = fit_to_data(int_times, mins[0])
  print("Scan 1 slope: ", popt[0]/1.e-6*100.e-12)

  popt1, pcov1 = fit_to_data(int_times, mins[1])
  print("Scan 2 slope: ", popt1[0]/1.e-6*100.e-12)

  popt2, pcov2 = fit_to_data(fit_times_10_50, np.asarray(mins[2])[fit_mask_10_50])
  print("Scan 3 slope: ", popt2[0]/1.e-6*100.e-12)

  popt3, pcov3 = fit_to_data(int_times, np.asarray(mins[3]))
  print("Scan 4 slope: ", popt3[0]/1.e-6*100.e-12)

  popt4, pcov4 = fit_to_data(fit_times_10_50, np.asarray(mins[4])[fit_mask_10_50])
  print("Scan 5 slope: ", popt4[0]/1.e-6*100.e-12)
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

  #axs1.plot(int_times*10., mins2[2], marker='o', label="2 ns cable")
  #axs1.plot(int_times*10., mins2[3], marker='o', label="4 ns cable")
  #axs1.plot(int_times*10., mins2[4], marker='o', label="8 ns cable")
  #axs1.plot(int_times*10., mins2[5], marker='o', label="16 ns cable")
  #axs1.plot(int_times*10., mins2[6], marker='o', label="32 ns cable")
  line, = axs1.plot(int_times, mins[0], marker='o', label="Scan 1: 1V, 100 kOhms")
  axs1.plot(int_times, linear_fit(int_times, popt[0], popt[1]), linestyle='--', color=line.get_color())

  line, = axs1.plot(int_times, mins[1], marker='o', label="Scan 2: 1V, 75 kOhms")
  axs1.plot(int_times, linear_fit(int_times, popt1[0], popt1[1]), linestyle='--', color=line.get_color())

  line, = axs1.plot(int_times, mins[2], marker='o', label="Scan 3: 1V, 50 kOhms")
  axs1.plot(fit_times_10_50, linear_fit(fit_times_10_50, popt2[0], popt2[1]), linestyle='--', color=line.get_color())

  line, = axs1.plot(int_times, mins[3], marker='o', label="Scan 4: 1.333V, 100 kOhms")
  axs1.plot(int_times, linear_fit(int_times, popt3[0], popt3[1]), linestyle='--', color=line.get_color())

  line, = axs1.plot(int_times, mins[4], marker='o', label="Scan 5: 2V, 100 kOhms")
  axs1.plot(fit_times_10_50, linear_fit(fit_times_10_50, popt4[0], popt4[1]), linestyle='--', color=line.get_color())

  axs1.grid()

  axs1.legend(fontsize=14)
  axs1.set_xlabel("Integration time (us)")
  axs1.set_ylabel("SWIC scanner signal (V)")
  fig1.suptitle("200 us wide pulse, different voltages and resistances")

  fig1.savefig("Figure_1.pdf")