
from matplotlib.lines import lineStyles
import matplotlib.pyplot as plt

# mean in y in cm
mean_x = [0.23201, 0.17573, 0.10745]
mean_y = [0.77513, 0.80061, 0.78321]

beamshift_sigma = [0,0.05,0.1]

plt.plot(beamshift_sigma,mean_x,marker='.',label='Mean in x')
plt.plot(beamshift_sigma,mean_y,marker='.',label='Mean in y')
plt.legend()
plt.xlabel('Proton beam shift in Y (cm)')
plt.ylabel('Histogram mean at MM1 (cm)')

plt.show()