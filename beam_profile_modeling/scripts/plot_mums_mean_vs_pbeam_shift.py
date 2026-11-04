"""
Plotting shifts in the x and y mean of histograms and KDEs at all 3 mums. 

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
# list corresponds to beamshift of 
#                   [nominal, 1 sigma, 2 sigma]
# taken from 100x96 'raw' corrected histograms from Mo
mm1_hist_x_100x96 = [0.2320, 0.1757, 0.1075]
mm1_hist_y_100x96 = [0.7751, 0.8006, 0.7832]
mm2_hist_x_100x96 = [0.2192, 0.0605, 0.1222]
mm2_hist_y_100x96 = [0.6218, 1.3220, 2.0420]
mm3_hist_x_100x96 = [0.2540, 0.0432, -0.0152]
mm3_hist_y_100x96 = [0.2965, 1.0920, 1.9570]

mm1_hist_x_80x80 = [0.1552,0.1327,0.0784]
mm1_hist_y_80x80 = [0.5123,0.5183,0.4307]
mm2_hist_x_80x80 = [0.0689,-0.0026,0.0992]
mm2_hist_y_80x80 = [0.3575,0.8871,1.3996]
mm3_hist_x_80x80 = [0.1609,0.0297,-0.0160]
mm3_hist_y_80x80 = [0.0762,0.7174,1.3000]

mm1_kde_x_80x80 = [0.1597,0.1353,0.0799]
mm1_kde_y_80x80 = [0.5182,0.5229,0.4347]
mm2_kde_x_80x80 = [0.0765,-0.0079,0.0963]
mm2_kde_y_80x80 = [0.3540,0.9050,1.3986]
mm3_kde_x_80x80 = [0.1655,0.0297,-0.0051]
mm3_kde_y_80x80 = [0.0959,0.7330,1.2949]


def full_hist_lin_and_quadratic_fit():
  """
  Plot the linear and quatratic fits of mums position
  in y vs proton beam pos in y to match the Yonehara plot.
  """

  m1, b1 = np.polyfit(pbeamshift_cm, mm1_hist_y_100x96, 1)
  m2, b2 = np.polyfit(pbeamshift_cm, mm2_hist_y_100x96, 1)
  m3, b3 = np.polyfit(pbeamshift_cm, mm3_hist_y_100x96, 1)

  coeffs1 = np.polyfit(pbeamshift_cm, mm1_hist_y_100x96, 2)
  coeffs2 = np.polyfit(pbeamshift_cm, mm2_hist_y_100x96, 2)
  coeffs3 = np.polyfit(pbeamshift_cm, mm3_hist_y_100x96, 2)
  model1 = np.poly1d(coeffs1)
  model2 = np.poly1d(coeffs2)
  model3 = np.poly1d(coeffs3)

  plt.scatter(pbeamshift_cm, mm1_hist_y_100x96, marker='o', label='MM1')
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

  plt.scatter(pbeamshift_cm, mm2_hist_y_100x96, marker='o', label='MM2')
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

  plt.scatter(pbeamshift_cm, mm3_hist_y_100x96, marker='o',label='MM3')
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
  print(f"Saved fig \n{FILE_NAME}")


def good_region_hist_and_kde():
  
  # histogram linear fits 
  m1h, b1h = np.polyfit(pbeamshift_cm, mm1_hist_y_80x80, 1)
  m2h, b2h = np.polyfit(pbeamshift_cm, mm2_hist_y_80x80, 1)
  m3h, b3h = np.polyfit(pbeamshift_cm, mm3_hist_y_80x80, 1)

  # kde linear fits
  m1k, b1k = np.polyfit(pbeamshift_cm, mm1_kde_y_80x80, 1)
  m2k, b2k = np.polyfit(pbeamshift_cm, mm2_kde_y_80x80, 1)
  m3k, b3k = np.polyfit(pbeamshift_cm, mm3_kde_y_80x80, 1)

  plt.scatter(pbeamshift_cm, mm1_hist_y_80x80, marker='o', color='royalblue', label=f"Hist MM1, slope={round(m1h,3)}")
  plt.plot(pbeamshift_smooth, 
    m1h * pbeamshift_smooth + b1h, 
    ls='-',
    color='royalblue',
    #label=f"Linear fit, slope={round(m1h,3)}"
  )

  plt.scatter(pbeamshift_cm, mm2_hist_y_80x80, marker='o', color='gold', label=f"Hist MM2, slope={round(m2h,3)}")
  plt.plot(pbeamshift_smooth,
    m2h * pbeamshift_smooth + b2h,
    ls='-',
    color='gold',
    #label=f"Linear fit, slope={round(m2h,3)}"
  )


  plt.scatter(pbeamshift_cm, mm3_hist_y_80x80, marker='o', color='mediumseagreen', label=f"Hist MM3, slope={round(m3h,3)}")
  plt.plot(pbeamshift_smooth,
    m3h * pbeamshift_smooth + b3h, 
    ls='-',
    color='mediumseagreen',
    #label=f"Linear fit, slope={round(m3h,3)}"
  )



  plt.scatter(pbeamshift_cm, mm1_kde_y_80x80, marker='.', color='blue', label=f"KDE MM1, slope={round(m1k,3)}")
  plt.plot(pbeamshift_smooth, 
    m1k * pbeamshift_smooth + b1k, 
    ls='--',
    color='blue',
    #label=f"Linear fit, slope={round(m1k,3)}"
  )

  plt.scatter(pbeamshift_cm, mm2_kde_y_80x80, marker='.', color='darkorange', label=f"KDE MM2, slope={round(m2k,3)}")
  plt.plot(pbeamshift_smooth,
    m2k * pbeamshift_smooth + b2k,
    ls='--',
    color='darkorange', 
    #label=f"Linear fit, slope={round(m2k,3)}"
  )

  plt.scatter(pbeamshift_cm, mm3_kde_y_80x80, marker='.', color='darkgreen', label=f"KDE MM3, slope={round(m3k,3)}")
  plt.plot(pbeamshift_smooth,
    m3k * pbeamshift_smooth + b3k, 
    ls='--',
    color='darkgreen',
    #label=f"Linear fit, slope={round(m3k,3)}"
  )

  handles, labels = plt.gca().get_legend_handles_labels()
  order = [0, 3, 1, 4, 2, 5]
  plt.legend([handles[i] for i in order], [labels[i] for i in order])
  plt.title("Mean in y (160x160 cm hists and KDEs) vs proton beam shift in y")
  plt.xlabel('Proton beam shift in y (cm)')
  plt.ylabel('Mean in y (cm)')

  plt.xticks(np.linspace(0,0.1,5))
  plt.ylim(-0.1,1.6)
  plt.yticks([0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50])

  plt.grid(True)
  plt.gca().set_axisbelow(True)


  print(f"MM1 slope difference: {np.abs(m1h-m1k):.4f}")
  print(f"MM2 slope difference: {np.abs(m2h-m2k):.4f}")
  print(f"MM3 slope difference: {np.abs(m3h-m3k):.4f}")

  # hist - kde mean difference at each proton beam shift, for each MM
  mean_diffs = {
    "MM1": np.array(mm1_hist_y_80x80) - np.array(mm1_kde_y_80x80),
    "MM2": np.array(mm2_hist_y_80x80) - np.array(mm2_kde_y_80x80),
    "MM3": np.array(mm3_hist_y_80x80) - np.array(mm3_kde_y_80x80),
  }

  max_diff = 0
  max_mm = None
  max_shift = None
  for mm, diffs in mean_diffs.items():
    for shift, diff in zip(pbeamshift_cm, diffs):
      print(f"{mm}, pbeam shift {shift} cm: hist - kde = {diff:.4f} cm")
      if np.abs(diff) > max_diff:
        max_diff = np.abs(diff)
        max_mm = mm
        max_shift = shift

  print(f"Max difference in mean: {max_diff:.4f} cm ({max_mm}, pbeam shift {max_shift} cm)")

  all_diffs = np.concatenate(list(mean_diffs.values()))
  print(f"Average difference in mean: {np.mean(all_diffs):.4f} cm")
  print(f"Average absolute difference in mean: {np.mean(np.abs(all_diffs)):.4f} cm")


  FILE_NAME = "mums_hist_and_kde_y_means_vs_pbeam_shift_pY.pdf"
  plt.savefig(os.path.join(OUTPUT_DIR, FILE_NAME))
  print(f"Saved fig \n{FILE_NAME}")

good_region_hist_and_kde()