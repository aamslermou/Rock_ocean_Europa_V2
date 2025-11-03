# This code computes the partial molar volume for a given gas at a given 
# temperature using interpolations formulas based on experimental data 
# For CO2 it is from : https://www.osti.gov/servlets/purl/790022 or Interpolation based on data from Rumpf et al 1993 a 
# For NH3 : Interpolation bsaed on data from Rumpf et al 1993 b 

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def partial_molar_volume(molecule, T):
    if molecule=='CO2':
        x = [313.15, 323.15, 333.15, 353.15, 393.15, 413.15, 433.15]
        y = [33.4, 34.0, 34.7, 36.3, 40.8, 43.8, 47.5]
        cs = CubicSpline(x, y)
        v = cs(T)
        # T = T -273
        # v = 37.51 - 9.585e-2*T + 8.740e-4*T**2 - 5.044e-7*T**3 
    elif molecule=='NH3':
        x = [298.15, 333.15, 353.15, 393.15, 413.15, 433.15]
        y = [30.0, 30.7, 32.1, 36.2, 38.9, 42.2]
        cs = CubicSpline(x, y)
        v = cs(T)
        # xs = np.linspace(298.15, 433.15)
        # fig,ax = plt.subplots()
        # ax.plot(x, y, 'o', label = 'data')
        # ax.plot(xs,cs(xs), label = 'extrapolation')
        # plt.show()
    return v*1e-6 

#print(partial_molar_volume('NH3', 340))

# T_span = np.linspace(273.15,360)

# plt.figure(1)

# plt.plot(T_span, partial_molar_volume('CO2',T_span))
# plt.plot(T_span, partial_molar_volume('NH3',T_span))

# plt.show()