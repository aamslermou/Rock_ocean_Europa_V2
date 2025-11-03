# This code displays de the functions used to compute the Henry constant for CO2/H2O and NH3/H2O couples.
# Due to the formulation of the liquid-vapor equilibrium, the henry constants have to be in the mole fraction scale  
# The formula used in this code are from Kawazuishi et al (1987)
# The formula used in this code are from Rumpf and Maurer (1993)
# The formula used for species other than H2O, CO2 and NH3 are from Warnek and Williams (2012)
# H [Pa]
# T [K]

import numpy as np
import matplotlib.pyplot as plt

Mw = 18.01528

def Henry_H2O_darde(T,component):
    if component == 'CO2':
        H = (1/Mw)*np.exp(192.876 - 9624.4/T + 1.441e-2*T - 28.749*np.log(T))
    elif component == 'NH3':
        H = (1/Mw)*np.exp(3.932 - 1879.02/T - 355134.1/(T**2))
    return H 

def Henry_H2O_bieling(T,component):
    if component == 'CO2':
        H = np.exp(-4.255 - 5678.0/T - 1.2191e6/(T**2))
    elif component == 'NH3':
        H = np.exp(3.932 - 1879.02/T - 355134.1/(T**2))
    return H 


def Henry_H2O_Kawasuishi_Warnek(T,component):
    if component == 'CO2':
        H = np.exp(-17060.71/T - 68.31596*np.log(T) + 0.06598907*T + 430.1920)
        # Conversion from [bar*kg/mol] to [Pa]
        H = H*1e5*1.01325*1e5/1845.40
    if component == 'NH3':
        H = np.exp(-7579.948/T - 13.58857*np.log(T) + 0.008596972*T + 96.23184)
        # Conversion from [bar*kg/mol] to [Pa]
        H = H*1e5*1.01325*1e5/1845.40
    if component == 'CH4':
        #H = 2.9477e6 - 44139*T + 246.83*T**2 - 0.64697*T**3 + 8.0669e-4*T**4 - 3.8747e-7
        H = np.exp(-211.28 + 10447.9/T + 29.780*np.log(T))
        # Conversion from [mol/(dm^3*atm)] to [Pa]
        H = 1.01325*1e5*55341.9/(H*1e3)
    if component == 'CO':
        H = np.exp(-178.0 + 8750/T + 24.875*np.log(T))
        # Conversion from [mol/(dm^3*atm)] to [Pa]
        #H = 1.01325*1e5/(1.83089*H*1e-3*1.01325*1e5)
        H = 1.01325*1e5*55341.9/(H*1e3)
    if component == 'N2':
        H = np.exp(-177.57 + 8632.1/T + 24.798*np.log(T))
        # Conversion from [mol/(dm^3*atm)] to [Pa]
        H = 1.01325*1e5*55341.9/(H*1e3)
    if component == 'Ar':
        H = np.exp(-146.40 + 7476.3/T + 20.140*np.log(T))
        # Conversion from [mol/(dm^3*atm)] to [Pa]
        H = 1.01325*1e5*55341.9/(H*1e3)
    if component == 'Kr':
        H = np.exp(-174.52 + 9101.7/T + 24.221*np.log(T))
        # Conversion from [mol/(dm^3*atm)] to [Pa]
        H = 1.01325*1e5*55341.9/(H*1e3)
    if component == 'Xe':
        H = np.exp(-197.21 + 10521/T + 27.466*np.log(T))
        # Conversion from [mol/(dm^3*atm)] to [Pa]
        H = 1.01325*1e5*55341.9/(H*1e3)
    if component == 'O2':
        H = np.exp(-173.33 + 8747.5/T + 24.453*np.log(T))
        # Conversion from [mol/(dm^3*atm)] to [Pa]
        H = 1.01325*1e5*55341.9/(H*1e3)
    return H 


def Henry_H2O_Rumpf(T,component):
    if component == 'CO2':
        H = np.exp(192.876 - 9624.4/T + 1.441e-2*T - 28.749*np.log(T))
        # Conversion from [Mpa*kg/mol] to [Pa]
        H = H*1e6*1.01325*1e5/1845.40
    elif component == 'NH3':
        H = np.exp(3.932 - 1879.02/T - 355134.1/(T**2))
        # Conversion from [Mpa*kg/mol] to [Pa]
        H = H*1e6*1.01325*1e5/1845.40
    return H 

# for T in np.linspace(275,386,50):
#     print(Henry_H2O_Kawasuishi_Warnek(T,'CH4'))

# for T in np.linspace(275,386,50):
# #     print(Henry_H2O_Kawasuishi_Warnek(T,'NH3'))

# T_span = np.linspace(273.15, 300, 50)

# plt.figure(1)

# plt.plot(T_span, Henry_H2O_Rumpf(T_span, 'NH3'))
# plt.plot(T_span, Henry_H2O_Rumpf(T_span, 'CO2'))

# plt.yscale('log')
# plt.show()

