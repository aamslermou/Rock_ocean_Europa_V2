# This code displays the method used to compute the constants of reactions K [mol/kg]
# They are computed using non linear relation from Kawazuishi et Prausnitz (1987) and Rumpf and Maurer (1993) for H2O :
# ln(K) = A1/T + A2*ln(T) +A3*T + A4   where Ai are constants from this paper and T is in K

import numpy as np

def K_computation_Kawazuishi(T,component):
        
    if component == 'H2O':
        #K = np.exp(-13445.9/T - 22.4773*np.log(T) - 0*T + 140) 
        K = 1e-14
    elif component == 'NH3':
        K = np.exp(-5914.082/T - 15.06399*np.log(T) - 0.01100801*T + 97.97152)
    elif component == 'CO2':
        K = np.exp(-7726.010/T - 14.50613*np.log(T) - 0.02798420*T + 102.2755)
    elif component == 'HCO3-':
        K = np.exp(-9137.258/T - 18.11192*np.log(T) - 0.02245619*T + 116.7371)
    else:
        K = np.exp(604.1164/T - 4.017263*np.log(T) +  0.005030950*T + 20.15214)

    return K

def K_computation_Edwards(T,component):
         
    if component == 'H2O':
        K = np.exp(-13445.9/T - 22.4773*np.log(T) - 0*T + 140) 
    elif component == 'NH3':
        K = np.exp(-3335.7/T + 1.4971*np.log(T) - 0.0370566*T + 2.76)
    elif component == 'CO2':
        K = np.exp(-12092.1/T  - 36.7816*np.log(T) - 0.0*T + 235.482)
    elif component == 'HCO3-':
        K = np.exp(-12431.7/T - 35.4819*np.log(T) - 0.0*T + 220.067)
    else:
        K = np.exp(-8.6 + 2900/T)

    return K

#print(K_computation_Kawazuishi(300,'HCO3-'))
#print(K_computation_Edwards(273,'H2O'),K_computation_Kawazuishi(298,'CO2'),K_computation_Kawazuishi(298,'NH3'),K_computation_Kawazuishi(298,'HCO3-'),K_computation_Kawazuishi(298,'g'))
#