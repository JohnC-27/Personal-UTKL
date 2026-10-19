"""
Plotting shifts in the x and y mean of histograms and KDEs at all 3 mums. 

Current state is plotting the linear and quatratic fits of mums position
in y vs proton beam pos in y to match the Yonehara plot.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")

# actual proton beam shift in cm [nominal, 1 sigma, 2 sigma]
pbeamshift_cm = [0,0.05,0.1]
# smooth for plotting fits
pbeamshift_smooth = np.linspace(pbeamshift_cm[0], pbeamshift_cm[-1], 200)

# shift in mean at each MM in cm, rounded to the nearest micron
# list corresponds to beamshift of [nominal, 1 sigma, 2 sigma]
mm1_hist_x = [0.2320, 0.1757, 0.1075]
mm1_hist_y = [0.7751, 0.8006, 0.7832]
mm2_hist_x = [0.2192, 0.0605, 0.1222]
mm2_hist_y = [0.6218, 1.3220, 2.0420]
mm3_hist_x = [0.2540, 0.0432, -0.0152]
mm3_hist_y = [0.2965, 1.0920, 1.9570]


m1, b1 = np.polyfit(pbeamshift_cm, mm1_hist_y, 1)
m2, b2 = np.polyfit(pbeamshift_cm, mm2_hist_y, 1)
m3, b3 = np.polyfit(pbeamshift_cm, mm3_hist_y, 1)

coeffs1 = np.polyfit(pbeamshift_cm, mm1_hist_y, 2)
coeffs2 = np.polyfit(pbeamshift_cm, mm2_hist_y, 2)
coeffs3 = np.polyfit(pbeamshift_cm, mm3_hist_y, 2)
model1 = np.poly1d(coeffs1)
model2 = np.poly1d(coeffs2)
model3 = np.poly1d(coeffs3)


plt.scatter(pbeamshift_cm, mm1_hist_y, marker='o', label='MM1')
plt.plot(pbeamshift_smooth, 
  m1 * pbeamshift_smooth + b1, 
  ls='-',
  #label=f"Linear fit, slope={round(m1,3)}"
)
plt.plot(pbeamshift_smooth, 
  model1(pbeamshift_smooth), 
  ls='--',
  color='tab:blue'
)

plt.scatter(pbeamshift_cm, mm2_hist_y, marker='o', label='MM2')
plt.plot(pbeamshift_smooth,
  m2 * pbeamshift_smooth + b2,
  ls='-', 
  #label=f"Linear fit, slope={round(m2,3)}"
)
plt.plot(pbeamshift_smooth, 
  model2(pbeamshift_smooth), 
  ls='--',
  color='tab:orange'
)

plt.scatter(pbeamshift_cm, mm3_hist_y, marker='o',label='MM3')
plt.plot(pbeamshift_smooth,
  m3 * pbeamshift_smooth + b3, 
  ls='-',
  #label=f"Linear fit, slope={round(m3,3)}"
)
plt.plot(pbeamshift_smooth, 
  model3(pbeamshift_smooth), 
  ls='--',
  color='tab:green'
)

plt.legend()
plt.title("Mean in y vs proton beam shift in y")
plt.xlabel('Proton beam shift in y (cm)')
plt.ylabel('Mean in y (cm)')

plt.xticks(np.linspace(0,0.1,5))
plt.ylim(-0.1,2.2)
plt.yticks(np.linspace(0,2,9))

plt.grid(True)
plt.gca().set_axisbelow(True)


FILE_NAME = "mums_y_hist_means_vs_pbeam_pY.pdf"
plt.savefig(os.path.join(OUTPUT_DIR, FILE_NAME))
print(f"Saved fig to \n{FILE_NAME}")