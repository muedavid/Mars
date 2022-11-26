import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from mars import Mars

emissivity_dessert = 0.5
emissivity_PV = 0.5
absorptivity_dessert = 0.5
absorptivity_PV = 0.5
delta_time = 24 * 60 * 60
f = 0.15
cp = 1
T_Atmosphere = 200
rho = 1
num_days = 700

# formula not given. Just a placeholder. replace values in mars.py
x = 1e-5

Mars = Mars(absorptivity_PV, absorptivity_dessert, emissivity_dessert, emissivity_PV, delta_time, f,
            cp, T_Atmosphere, rho, x)

# Just any simulated values. Replace by real ones. Numpy vector and each entry is average of the given day
L_in = np.cos(np.linspace(0, 2 * np.pi, num_days)) * 5 + 5
S_in = np.cos(np.linspace(0, 2 * np.pi, num_days)) * 5 + 5
r_H_dessert = np.ones(shape=num_days)
r_H_PV = r_H_dessert / 2

Temperature_init = np.ones(shape=(2 * num_days,)) * 273
root = fsolve(lambda Temperature: Mars.system(Temperature, num_days=num_days, l_in=L_in, s_in=S_in,
                                              r_H_PV=r_H_PV, r_H_dessert=r_H_dessert),
              Temperature_init)

plt.figure()
plt.plot(np.arange(0, num_days), root[0:num_days])
plt.plot(np.arange(0, num_days), root[num_days:2 * num_days])
plt.show()
