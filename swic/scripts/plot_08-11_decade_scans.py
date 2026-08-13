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
  data_dir = os.path.join(os.path.dirname(__file__), "..", "2026-08-11")
  
  filelist1 = []
  for files in os.listdir(data_dir):
    if files.endswith(".csv"):
      filelist1.append(os.path.join(data_dir, files))
  
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
  
  #print(filelist1)
  #print(split_channels_numpy(filelist[0]))
  mins_100 = [get_mins(split_channels_numpy(file))[1] for file in filelist] #[np.min(channels['set_1'],axis=1), np.min(channels['set_2'],axis=1), np.min(channels['set_3'],axis=1), np.min(channels['set_4'],axis=1)]
  maxes = [get_max(split_channels_numpy(file))[1] for file in filelist]

  
  #mins2 = [get_mins(split_channels_numpy(file))[1] for file in filelist1]

  mins_1000 = [get_mins(split_channels_numpy(file))[3] for file in filelist]

  mins_10000 = [get_mins(split_channels_numpy(file))[2] for file in filelist]

  mins_100000 = [get_mins(split_channels_numpy(file))[0] for file in filelist]

  popt, pcov = fit_to_data(int_times_2, mins_100[1])
  print(popt[0]/1.e-6*100.e-12)

  popt1, pcov1 = fit_to_data(int_times_2, mins_100[2])
  print(popt1[0]/1.e-6*100.e-12)

  popt2, pcov2 = fit_to_data(int_times_2[1:], mins_100[5][1:])
  print(popt2[0]/1.e-6*100.e-12)

  popt3, pcov3 = fit_to_data(int_times_2[1:], mins_100[5][1:])
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

  linstyles_for_current = ['-','--',':']

  labels = ["10 uA", "20 uA", "50 uA", "100 uA"]

  for i in range(8):
    if i<4:
      axs1.plot(int_times*10., mins_100[i], marker='o', color=current_colors_tp[i%4],linestyle=linstyles_for_current[0],label=labels[i])
    else:
      axs1.plot(int_times*10., mins_100[i], marker='o', color=current_colors_tp[i%4],linestyle=linstyles_for_current[1])

  axs1.grid()

  axs1.legend(fontsize=14)
  axs1.set_xlabel("Integration time (us)")
  axs1.set_ylabel("SWIC scanner signal (V)")
  axs1.set_ylim(-10.5,0)
  fig1.suptitle("30 us wide pulses into 100 pF board")

  fig1.savefig(os.path.join(OUTPUT_DIR, "100pF_board_08-11.pdf"))

  fig2, axs2 = plt.subplots(figsize=(10,8))

  labels = ["100 uA", "200 uA", "500 uA", "1 mA"]

  for i in range(8,16):
    if i<12:
      axs2.plot(int_times*10., mins_1000[i], marker='o', color=current_colors_tp[i%4],linestyle=linstyles_for_current[0],label=labels[i%4])
    else:
      axs2.plot(int_times*10., mins_1000[i], marker='o', color=current_colors_tp[i%4],linestyle=linstyles_for_current[1])
    #axs2.plot(int_times*10., mins_1000[i], marker='o', color=current_colors_tp[i%4],linestyle=linstyles_for_current[2-((i%12)//4)])

  #axs2.plot(int_times*10., mins_1000[37], marker='o', color='black')

  axs2.grid()

  axs2.legend(fontsize=14)
  axs2.set_xlabel("Integration time (us)")
  axs2.set_ylabel("SWIC scanner signal (V)")
  axs2.set_ylim(-8.5,0)
  fig2.suptitle("30 us wide pulses into 1000 pF board")

  fig2.savefig(os.path.join(OUTPUT_DIR, "1000pF_board_08-11.pdf"))

  fig3, axs3 = plt.subplots(figsize=(10,8))

  labels = ["1 mA", "2 mA", "5 mA", "10 mA"]

  for i in range(16,24):
    #if i==35 or i==31 or i==27:
    #  continue
    if i<20:
      axs3.plot(int_times*10., mins_10000[i], marker='o', color=current_colors_tp[i%4],linestyle=linstyles_for_current[0],label=labels[i%4])
    else:
      axs3.plot(int_times*10., mins_10000[i], marker='o', color=current_colors_tp[i%4],linestyle=linstyles_for_current[1])

  #axs3.plot(int_times*10., mins_10000[38], marker='o', color=current_colors_tp[3], linestyle=linstyles_for_current[0])
  #axs3.plot(int_times*10., mins_10000[39], marker='o', color=current_colors_tp[3], linestyle=linstyles_for_current[1])
  #axs3.plot(int_times*10., mins_10000[40], marker='o', color=current_colors_tp[3], linestyle=linstyles_for_current[2])

  #axs3.plot(int_times*10., mins_10000[36], marker='o', color='black',label="1 mA, different resistor")

  axs3.grid()

  axs3.legend(fontsize=14)
  axs3.set_xlabel("Integration time (us)")
  axs3.set_ylabel("SWIC scanner signal (V)")
  axs3.set_ylim(-3.5,0)
  fig3.suptitle("30 us wide pulses into 10000 pF board")

  fig3.savefig(os.path.join(OUTPUT_DIR, "10000pF_board_08-11.pdf"))

  labels = ["1 mA", "2 mA", "5 mA", "10 mA"]

  fig4, axs4 = plt.subplots(figsize=(10,8))
  axs4.plot(int_times*10., mins_100000[24], marker='o', color=current_colors_tp[0],linestyle=linstyles_for_current[0],label="10 mA")
  axs4.plot(int_times*10., mins_100000[25], marker='o', color=current_colors_tp[1],linestyle=linstyles_for_current[0],label="20 mA")
  axs4.plot(int_times*10., mins_100000[26], marker='o', color=current_colors_tp[2],linestyle=linstyles_for_current[0],label="50 mA")
  axs4.plot(int_times*10., mins_100000[27], marker='o', color=current_colors_tp[0],linestyle=linstyles_for_current[1])
  axs4.plot(int_times*10., mins_100000[28], marker='o', color=current_colors_tp[1],linestyle=linstyles_for_current[1])
  axs4.plot(int_times*10., mins_100000[29], marker='o', color=current_colors_tp[2],linestyle=linstyles_for_current[1])
  #axs4.plot(int_times*10., mins_100000[30], marker='o', color=current_colors_tp[0],linestyle=linstyles_for_current[2])
  #axs4.plot(int_times*10., mins_100000[48], marker='o', color=current_colors_tp[1],linestyle=linstyles_for_current[2])
  #axs4.plot(int_times*10., mins_100000[49], marker='o', color=current_colors_tp[2],linestyle=linstyles_for_current[2])

  axs4.grid()

  axs4.legend(fontsize=14)
  axs4.set_xlabel("Integration time (us)")
  axs4.set_ylabel("SWIC scanner signal (V)")
  axs4.set_ylim(-0.5,0)
  fig4.suptitle("30 us wide pulses into 100000 pF board")

  fig4.savefig(os.path.join(OUTPUT_DIR, "100000pF_board_08-11.pdf"))