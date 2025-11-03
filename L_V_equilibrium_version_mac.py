# %%                                      ############## Functions importations ################

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from labellines import labelLine, labelLines

from gamma import gamma_UNIQUAC_2species
from gamma import gamma_UNIQUAC_Vext
from gamma import fct_gamma_CH4

from Henry import Henry_H2O_Kawasuishi_Warnek
from Henry import Henry_H2O_Rumpf
from Constants_reactions import K_computation_Kawazuishi 
from fugacity import FUG_COEFF_PRG_EOS
from partial_molar_volumes import partial_molar_volume
from Saturation_pressure import P_sat_water

from magnitude import magnitude

import os
import importlib

import Compo_ocean_Moon as comp
comp = importlib.reload(comp)


M_H2O     = 18.01528e-3 # [kg/mol]
rho_water = 1000        # [kg/m3]
R         = 8.314       # [J/(mol*K)]

##                                  ############### UNIQUAC importation ################

# import ctypes

# # Order of molecules : H2O, NH3, CO2, O2, N2, NH4+, H+, OH-, CO32-, HCO3-, NH2COO-

# #Define the input variables
# TKEL = ctypes.c_double(273.15)
# FEED = (ctypes.c_double * 23)(55.5, 0.966, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# # Define the other input variables with appropriate ctypes data types
# PBUB = ctypes.c_double(0)
# PRES = (ctypes.c_double * 50)()
# NGAS = ctypes.c_long(0)
# NBR = (ctypes.c_long * 50)()
# HGABS = (ctypes.c_double * 50)()
# AW = ctypes.c_double(0)
# GM = ctypes.c_double(0)
# ERROR = ctypes.c_long(0)

# #Load the DLL
# aqsol_dll = ctypes.windll.LoadLibrary(r'C:\Users\aamsler\Desktop\Uniquac\AQSOL25264.DLL')

# #Define the function prototype using the header information
# BUBBLEP = aqsol_dll.BUBBLEP
# BUBBLEP.argtypes = [
#     ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
#     ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),
#     ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
#     ctypes.POINTER(ctypes.c_long)
# ]

# BUBBLEP.restype = None
# %%
def L_V_equilibrium(T, nameSpecies, M_bulk, initial_guess, composition_liquidphase, partial_pressures, i, d_Ocean, R_Europa, g_Europa, u0_coeff, ut_coeff, k_coeff, species_phi, species):
    """
    Function to compute the equilibrium between the liquid and vapor phases of a system.
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
    nameSpecies : list
        List of species names.
    M_bulk : DataFrame
        Array containing the bulk composition of the hydrosphere.
    initial_guess : list
        Initial guess for the partial pressures of the species.
    composition_liquidphase : dataFrame
        Array to store the composition of the liquid phase.
    partial_pressures : dataFrame
        Array to store the partial pressures of the species.
    i : int
        Index for the current iteration.
    d_Ocean : float
        Depth of the ocean in meters.
    R_Europa : float
        Radius of Europa in meters.
    g_Europa : float
        Gravitational acceleration on Europa in m/s^2.
    u0_coeff : DataFrame
        Coefficients for the UNIQUAC model.
    ut_coeff : DataFrame
        Coefficients for the UNIQUAC model. 
    k_coeff : DataFrame
        Coefficients for the fugacity model.
    species_phi : dataFrame
        Species for the fugacity model.
    species : dataFrale
        Species constants.
    Returns
    -------
    composition_liquidphase, partial_pressures at iteration i
    """

    NN = 9
    psi = np.zeros((NN,NN))
    u   = np.zeros((NN,NN))
        
    for k in range(0,NN):
        for j in range(0,NN):
            u[k,j]   = u0_coeff.iat[k,j] + ut_coeff.iat[k,j]*(T-298.15)   # u computed with data from litterature (Darde et al. 2010)
        
    for k in range(NN):
        for j in range(NN):      
            psi[k,j] = np.exp(-(u[k,j]-u[j,j])/T)                         # Computation of psi from UNIQUAC model

    # Saturation coefficient of water : 
    PsatH2O = P_sat_water(T)   # [Pa]    

    VOcean = (4/3)*np.pi*(R_Europa**3 - (R_Europa - d_Ocean)**3)
    # M_bulk = np.array([comp.compo_67P[mol]*(VOcean*rho_water) for mol in nameSpecies]) # [kg]
    # print(M_bulk)
#          # Henry's constants
    H_CO2_H2O = Henry_H2O_Rumpf(T,'CO2')
    H_NH3_H2O = Henry_H2O_Rumpf(T,'NH3')
    H_CH4_H2O = Henry_H2O_Kawasuishi_Warnek(T,'CH4')
    H_Ar_H2O  = Henry_H2O_Kawasuishi_Warnek(T,'Ar')
    H_Kr_H2O  = Henry_H2O_Kawasuishi_Warnek(T,'Kr')
    H_Xe_H2O  = Henry_H2O_Kawasuishi_Warnek(T,'Xe')

        # Constants of reactions : 
    K_H2O     = K_computation_Kawazuishi(T,'H2O')
    K_CO2     = K_computation_Kawazuishi(T,'CO2')
    K_NH3     = K_computation_Kawazuishi(T,'NH3')
    K_HCO3m   = K_computation_Kawazuishi(T,'HCO3-')
    K_NH2COOm = K_computation_Kawazuishi(T,'NH2COO-')

    H = [H_CO2_H2O, H_NH3_H2O, H_CH4_H2O, H_Ar_H2O, H_Kr_H2O, H_Xe_H2O]

        # Partial molar volumes at infinite dilution in water
    if T < 313.15:
        v_CO2 = 0
        v_NH3 = 0
    elif 313.15 < T < 333.15:
        v_CO2 = partial_molar_volume('CO2',T)
        v_NH3 = 0    
    else:
        v_CO2 = partial_molar_volume('CO2',T)
        v_NH3 = partial_molar_volume('NH3',T)

               ############# Implementation of the system's functions ################

## Definition of the anonymous functions to be solved : 

    # Equality between fugacities in liquid and gaseous phases : 

    #f1_H2O = lambda P,y_H2O,x_H2O,phi_H2O,gamma_H2O,v : (phi_H2O*P*y_H2O)/(gamma_H2O*x_H2O*PsatH2O*np.exp((v*(P-PsatH2O))/(R*T))) - 1
    f1_H2O_y = lambda P,y_H2O,x_H2O,phi_H2O,gamma_H2O   : (gamma_H2O*x_H2O*PsatH2O)/(phi_H2O*P*y_H2O) - 1
    f2_CO2_y = lambda P,y_CO2,x_CO2,phi_CO2,gamma_CO2   : (gamma_CO2*x_CO2*H_CO2_H2O*np.exp((v_CO2*(P-PsatH2O))/(R*T)))/(phi_CO2*P*y_CO2) - 1
    f3_NH3_y = lambda P,y_NH3,x_NH3,phi_NH3,gamma_NH3   : (gamma_NH3*x_NH3*H_NH3_H2O*np.exp((v_NH3*(P-PsatH2O))/(R*T)))/(phi_NH3*P*y_NH3) - 1  
    f5_CH4_y = lambda P,y_CH4,x_CH4,phi_CH4,gamma_CH4   : (gamma_CH4*x_CH4*H_CH4_H2O)/(phi_CH4*P*y_CH4) - 1
    f7_Ar_y  = lambda P,y_Ar,x_Ar,phi_Ar                : (x_Ar*H_Ar_H2O)/(phi_Ar*P*y_Ar) - 1
    f8_Kr_y  = lambda P,y_Kr,x_Kr,phi_Kr                : (x_Kr*H_Kr_H2O)/(phi_Kr*P*y_Kr) - 1
    f9_Xe_y  = lambda P,y_Xe,x_Xe,phi_Xe                : (x_Xe*H_Xe_H2O)/(phi_Xe*P*y_Xe) - 1

    f1_H2O_x = lambda P,y_H2O,x_H2O,phi_H2O,gamma_H2O   : (phi_H2O*P*y_H2O)/(gamma_H2O*x_H2O*PsatH2O) - 1
    f2_CO2_x = lambda P,y_CO2,x_CO2,phi_CO2,gamma_CO2   : (phi_CO2*P*y_CO2)/(gamma_CO2*x_CO2*H_CO2_H2O*np.exp((v_CO2*(P-PsatH2O))/(R*T))) - 1
    f3_NH3_x = lambda P,y_NH3,x_NH3,phi_NH3,gamma_NH3   : (phi_NH3*P*y_NH3)/(gamma_NH3*x_NH3*H_NH3_H2O*np.exp((v_NH3*(P-PsatH2O))/(R*T))) - 1  
    f5_CH4_x = lambda P,y_CH4,x_CH4,phi_CH4,gamma_CH4   : (phi_CH4*P*y_CH4)/(gamma_CH4*x_CH4*H_CH4_H2O) - 1
    f7_Ar_x  = lambda P,y_Ar,x_Ar,phi_Ar                : (phi_Ar*P*y_Ar)/(x_Ar*H_Ar_H2O) - 1
    f8_Kr_x  = lambda P,y_Kr,x_Kr,phi_Kr                : (phi_Kr*P*y_Kr)/(x_Kr*H_Kr_H2O) - 1
    f9_Xe_x  = lambda P,y_Xe,x_Xe,phi_Xe                : (phi_Xe*P*y_Xe)/(x_Xe*H_Xe_H2O) - 1   

    f1_H2O = lambda P_H2O,x_H2O,phi_H2O,gamma_H2O     : (phi_H2O*P_H2O)/(gamma_H2O*x_H2O*PsatH2O) - 1
    f2_CO2 = lambda P_CO2,x_CO2,phi_CO2,gamma_CO2,Ptot   : (phi_CO2*P_CO2)/(gamma_CO2*x_CO2*H_CO2_H2O*np.exp((v_CO2*(Ptot-PsatH2O))/(R*T))) - 1
    f3_NH3 = lambda P_NH3,x_NH3,phi_NH3,gamma_NH3,Ptot   : (phi_NH3*P_NH3)/(gamma_NH3*x_NH3*H_NH3_H2O*np.exp((v_NH3*(Ptot-PsatH2O))/(R*T))) - 1  

    # Total of molar fractions : 

    f4      = lambda x_H2O, x_CO2, x_NH3, x_CH4, x_Ar, x_Kr, x_Xe : x_H2O + x_CO2 + x_NH3 + x_CH4 + x_Ar + x_Kr + x_Xe - 1
    f4_bis  = lambda x_H2O,x_CO2,x_NH3,x_HCO3m,x_Hp,x_NH4p,x_CO32m,x_OHm,x_NH2COOm, x_CH4, x_Ar, x_Kr, x_Xe : x_H2O + x_CO2 + x_NH3 + x_HCO3m + x_Hp +\
            x_NH4p + x_CO32m + x_OHm + x_NH2COOm + x_CH4  + x_Ar + x_Kr + x_Xe - 1
    f5      = lambda y_H2O, y_CO2, y_NH3, y_CH4, y_Ar, y_Kr, y_Xe : y_H2O + y_CO2 + y_NH3 + y_CH4 + y_Ar + y_Kr + y_Xe - 1

    f14    = lambda x_H2O, m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm : (1000/M_H2O)/(1000/M_H2O + m_CO2 + m_NH3 + m_HCO3m + m_Hp + m_NH4p + m_CO32m + m_OHm + m_NH2COOm) * 1/x_H2O - 1

    # Mass balance : 

    f6     = lambda m_CO2,m_HCO3m,m_NH2COOm,m_CO32m  : (m_CO2 + m_HCO3m + m_NH2COOm + m_CO32m)/m_CO2_tot - 1
    f7     = lambda m_NH3,m_NH4p, m_NH2COOm          : (m_NH3 + m_NH4p + m_NH2COOm)/m_NH3_tot - 1 

    # Charge balance :

    f8     = lambda m_HCO3m,m_Hp,m_NH4p,m_CO32m,m_OHm,m_NH2COOm : (m_Hp + m_NH4p)/(m_HCO3m + 2*m_CO32m + m_OHm + m_NH2COOm) - 1

    # Constants of reactions equations : 

    f9     = lambda x_H2O,m_OHm,m_Hp,gamma_Hp,gamma_OHm,gamma_H2O                               : (m_OHm*m_Hp)/(x_H2O) * (gamma_OHm*gamma_Hp/(gamma_H2O)) * 1/K_H2O - 1

    f10    = lambda x_H2O,m_CO2,m_HCO3m,m_Hp,gamma_H2O,gamma_CO2,gamma_HCO3m,gamma_Hp           : (m_HCO3m*m_Hp/(m_CO2*x_H2O)) * (gamma_HCO3m*gamma_Hp/(gamma_CO2*gamma_H2O)) * 1/K_CO2 - 1

    f11    = lambda x_H2O,m_HCO3m,m_Hp,m_CO32m,gamma_HCO3m,gamma_Hp,gamma_CO32m                 : (m_CO32m*m_Hp/(m_HCO3m*x_H2O)) * (gamma_CO32m*gamma_Hp/(gamma_HCO3m*x_H2O)) * 1/K_HCO3m - 1

    f12    = lambda x_H2O,m_NH3,m_NH4p,m_OHm,gamma_H2O,gamma_NH3,gamma_NH4p,gamma_OHm           : (m_NH4p*m_OHm/(m_NH3*x_H2O)) * (gamma_NH4p*gamma_OHm/(gamma_NH3*gamma_H2O)) * 1/K_NH3 - 1

    f13    = lambda x_H2O,m_NH3,m_HCO3m,m_NH2COOm,gamma_H2O,gamma_NH3,gamma_HCO3m,gamma_NH2COOm : (m_NH2COOm*x_H2O/(m_NH3*m_HCO3m)) * (gamma_NH2COOm*gamma_H2O/(gamma_NH3*gamma_HCO3m)) * 1/K_NH2COOm - 1

##
    def myfun1(z):
        
        x_H2O, x_CO2, x_NH3, x_CH4, x_Ar, x_Kr, x_Xe  = z

        # Fugacities in gaseous phase : 

        y   = [y_H2O, y_CO2, y_NH3, y_CH4, y_Ar, y_Kr, y_Xe]

        phi = FUG_COEFF_PRG_EOS(R,T,P,species_phi,k_coeff,y)[0]

        phi_H2O, phi_CO2, phi_NH3, phi_CH4, phi_Ar, phi_Kr, phi_Xe = phi

        # Coefficients of activity : 

        gamma_H2O, gamma_CO2, gamma_NH3 = [1,1,1]
        gamma_CH4 = fct_gamma_CH4(T,x_CH4)

        return([f2_CO2_y(P, y_CO2, x_CO2, phi_CO2, gamma_CO2), f3_NH3_y(P, y_NH3, x_NH3, phi_NH3, gamma_NH3),\
        f5_CH4_y(P, y_CH4, x_CH4, phi_CH4, gamma_CH4), f7_Ar_y(P, y_Ar, x_Ar, phi_Ar), f8_Kr_y(P, y_Kr, x_Kr, phi_Kr),\
            f9_Xe_y(P, y_Xe, x_Xe, phi_Xe), f4(x_H2O, x_CO2, x_NH3, x_CH4, x_Ar, x_Kr, x_Xe)])


    def Chemical_eq_fun(z):

        x_H2O, m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm = z

        # Fugacities in gaseous phase : 

        # y = [y_H2O, y_CO2, y_NH3, y_O2, y_CO, y_CH4, y_N2, y_Ar, y_Kr, y_Xe]

        # phi,Z = FUG_COEFF_PRG_EOS(R,T,P,species_phi,k_coeff,y)

        # phi_H2O, phi_CO2, phi_NH3, phi_O2, phi_CO, phi_CH4, phi_N2, phi_Ar, phi_Kr, phi_Xe = phi

        #print(phi)
        # Coefficient of activity : 
        x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm = np.array([m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm])*M_H2O*x_H2O

        gamma = gamma_UNIQUAC_Vext(T,x_H2O,x_CO2,x_NH3,x_HCO3m,x_Hp,x_NH4p,x_CO32m,x_OHm,x_NH2COOm,species,u0_coeff,ut_coeff)

        #print(gamma)
        gamma_H2O, gamma_CO2, gamma_NH3, gamma_HCO3m, gamma_Hp, gamma_NH4p, gamma_CO32m, gamma_OHm, gamma_NH2COOm = gamma

        #gamma_CH4 = fct_gamma_CH4(T,x_CH4)
        #gamma_CH4 = 1

        #print(gamma_CH4)
        return([f9(x_H2O,m_OHm,m_Hp,gamma_Hp,gamma_OHm,gamma_H2O),f10(x_H2O,m_CO2,m_HCO3m,m_Hp,gamma_H2O,gamma_CO2,gamma_HCO3m,gamma_Hp),f11(x_H2O,m_HCO3m,m_Hp,m_CO32m,gamma_HCO3m,gamma_Hp,gamma_CO32m),\
                f12(x_H2O,m_NH3,m_NH4p,m_OHm,gamma_H2O,gamma_NH3,gamma_NH4p,gamma_OHm),f13(x_H2O,m_NH3,m_HCO3m,m_NH2COOm,gamma_H2O,gamma_NH3,gamma_HCO3m,gamma_NH2COOm),\
                f14(x_H2O, m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm), f6(m_CO2,m_HCO3m,m_NH2COOm,m_CO32m), f7(m_NH3,m_NH4p, m_NH2COOm), f8(m_HCO3m,m_Hp,m_NH4p,m_CO32m,m_OHm,m_NH2COOm)])

    def L_V_eq_fun(z):

        P_H2O, P_CO2, P_NH3 = z

        # Fugacities in gaseous phase :

        Ptot = P_H2O + P_CO2 + P_NH3
        y = [P_H2O/Ptot, P_CO2/Ptot, P_NH3/Ptot]

        phi,Z = FUG_COEFF_PRG_EOS(R,T,Ptot,species_phi,k_coeff,y)

        phi_H2O, phi_CO2, phi_NH3 = phi
        #print(phi)
            # Coefficient of activity : 
        x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm = np.array([m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm])*M_H2O*x_H2O

        gamma = gamma_UNIQUAC_Vext(T,x_H2O,x_CO2,x_NH3,x_HCO3m,x_Hp,x_NH4p,x_CO32m,x_OHm,x_NH2COOm,species,u0_coeff,ut_coeff)

        #print(gamma)
        gamma_H2O, gamma_CO2, gamma_NH3, gamma_HCO3m, gamma_Hp, gamma_NH4p, gamma_CO32m, gamma_OHm, gamma_NH2COOm = gamma
        gamma_H2O, gamma_CO2, gamma_NH3 = [1,1,1]
        return([f1_H2O(P_H2O,x_H2O,phi_H2O,gamma_H2O),f2_CO2(P_CO2,x_CO2,phi_CO2,gamma_CO2,Ptot),f3_NH3(P_NH3,x_NH3,phi_NH3,gamma_NH3,Ptot)])

    def myfun3(z):
        
        x_H2O, x_CO2, x_CH4, x_Ar, x_Kr, x_Xe  = z

        # Fugacities in gaseous phase : 

        y   = [y_H2O, y_CO2, y_CH4, y_Ar, y_Kr, y_Xe]

        phi = FUG_COEFF_PRG_EOS(R,T,P,species,k_coeff,y)[0]

        phi_H2O, phi_CO2, phi_CH4, phi_Ar, phi_Kr, phi_Xe = phi

        # Coefficients of activity : 

        gamma_H2O, gamma_CO2 = [1,1]
        gamma_CH4 = fct_gamma_CH4(T,x_CH4)

        return([f2_CO2_y(P, y_CO2, x_CO2, phi_CO2, gamma_CO2),\
        f5_CH4_y(P, y_CH4, x_CH4, phi_CH4, gamma_CH4), f7_Ar_y(P, y_Ar, x_Ar, phi_Ar), f8_Kr_y(P, y_Kr, x_Kr, phi_Kr),\
            f9_Xe_y(P, y_Xe, x_Xe, phi_Xe), f4(x_H2O, x_CO2, x_CH4, x_Ar, x_Kr, x_Xe)])


                                ### Main program ###
##
    epsilon_high = np.zeros(6, dtype=np.float32)
    epsilon_low  = np.zeros(6, dtype=np.float32)

    for j, M in enumerate(M_bulk[1:]):
        if i != 2:
            if M < 1e13:
                epsilon_high[j] = 1e-3
                epsilon_low[j]  = 5e-4
            elif 1e13 < M < 1e15:
                epsilon_high[j] = 1e-2
                epsilon_low[j]  = 4e-3
            elif 1e15 < M < 5e16:
                epsilon_high[j] = 5
                epsilon_low[j]  = 1e-1
            elif 5e16 < M < 1e17:
                epsilon_high[j] = 100
                epsilon_low[j]  = 1
            elif 1e17 < M < 1e21:
                epsilon_high[j] = 1e3
                epsilon_low[j]  = 1e2
            elif 1e21 < M < 1e25:
                epsilon_high[j] = 1e4
                epsilon_low[j]  = 1e3
        else:
            if M < 1e17:
                epsilon_high[j] = 1e-2
                epsilon_low[j]  = 1e-3
            elif 1e17 < M :
                epsilon_high[j] = 1e0
                epsilon_low[j]  = 1e-2

   # epsilon = np.array([100, 100, 100, 10, 5, 1], dtype=np.float32)

        partial_pressures['H2O'][i] = P_sat_water(T)
        partial_pressures['CO2'][i] = initial_guess[0]
        partial_pressures['NH3'][i] = initial_guess[1]
        partial_pressures['CH4'][i] = initial_guess[2]
        partial_pressures['Ar'][i]  = initial_guess[3]
        partial_pressures['Kr'][i]  = initial_guess[4]
        partial_pressures['Xe'][i]  = initial_guess[5]

    all_met = True

    while True:
        # Calculate total pressure
        P = (partial_pressures['H2O'][i] + 
            partial_pressures['CO2'][i] + 
            partial_pressures['NH3'][i] +
            partial_pressures['CH4'][i] +
            partial_pressures['Ar'][i]  + 
            partial_pressures['Kr'][i]  +
            partial_pressures['Xe'][i])
        
        # Calculate molar fractions
        y_comp = np.array([partial_pressures[mol][i] / P for mol in ['H2O', 'CO2', 'NH3', 'CH4', 'Ar', 'Kr', 'Xe']])
        y_H2O, y_CO2, y_NH3, y_CH4, y_Ar, y_Kr, y_Xe = y_comp

        x0 = []
        NN = len(y_comp)
        
        for j in range(NN):
            if j == 0:
                    x0.append(0.99)
            else:
                    order_magnitude = magnitude(( P*y_comp[j]) / H[j-1])
                    x0.append(10 ** order_magnitude)

        # Solve using least squares
        sol = least_squares(myfun1, x0, method='lm', ftol=1e-10, xtol=1e-10)

        # Unpack solution
        x_H2O, x_CO2, x_NH3, x_CH4, x_Ar, x_Kr, x_Xe = sol.x
        res = sol.x

        # Calculate atmospheric and ocean masses
        M_atm   = np.array([(4 * np.pi * R_Europa**2 * partial_pressures[mol][i]) / g_Europa for mol in ['CO2', 'NH3', 'CH4', 'Ar', 'Kr', 'Xe']])
        M_ocean = np.array([(VOcean * 1000) * x for x in res[1:]])       # VOcean is in m^3

        # Update tolerances
        tolerances = (M_atm + M_ocean) - M_bulk[1:]

        # Flag to check if all tolerances are met
        all_met = True
        
        # Update partial pressures based on tolerances
        for k, mol in enumerate(['CO2', 'NH3', 'CH4', 'Ar', 'Kr', 'Xe']):

            if mol == 'CO2'  :                      # For CO2

                if np.abs(tolerances[k]) > 1e-3 * M_bulk[k]:  # Check if tolerance is still unmet
                    print(f"Checking {mol}: Tolerance = {tolerances[k]}, Condition = {tolerances[k] > 0.1 * M_bulk[k]}")
                    all_met = False  # Set flag to false if any tolerance is not met              
                        # Update epsilon based on the tolerance condition
                    epsilon = epsilon_low[k] if np.abs(tolerances[k]) < 0.01 * M_bulk[k] else epsilon_high[k]
                        
                        # Update partial pressures based on tolerance, skip if already met
                    if tolerances[k] > 0:
                        partial_pressures.at[i, mol] += (-1) * epsilon  # Decrease guessed pressure
                    else:
                        partial_pressures.at[i, mol] += epsilon         # Increase guessed pressure

            elif mol == 'CH4':       # For CH4

                if np.abs(tolerances[k]) > 1e-3 * M_bulk[k]:
                    print(f"Checking {mol}: Tolerance = {tolerances[k]}, Condition = {tolerances[k] > 0.1 * M_bulk[k]}")
                    all_met = False
                    epsilon = epsilon_low[k] if np.abs(tolerances[k]) < 0.01 * M_bulk[k] else epsilon_high[k]
                    if tolerances[k] > 0:
                        partial_pressures.at[i, mol] += (-1) * epsilon
                    else:
                        partial_pressures.at[i, mol] += epsilon

            elif mol == 'NH3':       # For NH3

                if np.abs(tolerances[k]) > 1e-3 * M_bulk[k]:
                    print(f"Checking {mol}: Tolerance = {tolerances[k]}, Condition = {tolerances[k] > 0.1 * M_bulk[k]}")
                    all_met = False
                    epsilon = epsilon_low[k] if np.abs(tolerances[k]) < 0.05 * M_bulk[k] else epsilon_high[k]
                    if M_bulk[2]> 1e16:
                        epsilon = epsilon/10
                    if tolerances[k] > 0:
                        partial_pressures.at[i, mol] += (-1) * epsilon
                    else:
                        partial_pressures.at[i, mol] += epsilon

            elif mol == 'Ar':       # For noble gases

                if np.abs(tolerances[k]) > 1e-3 * M_bulk[k]:
                    print(f"Checking {mol}: Tolerance = {tolerances[k]}, Condition = {tolerances[k] > 0.1 * M_bulk[k]}")
                    all_met = False
                    epsilon = epsilon_low[k] if np.abs(tolerances[k]) < 0.01 * M_bulk[k] else epsilon_high[k]
                    if tolerances[k] > 0:
                        partial_pressures.at[i, mol] += (-1) * epsilon
                    else:
                        partial_pressures.at[i, mol] += epsilon

            elif mol == 'Kr':       # For noble gases     

                if np.abs(tolerances[k]) > 1e-3 * M_bulk[k]:
                    print(f"Checking {mol}: Tolerance = {tolerances[k]}, Condition = {tolerances[k] > 0.1 * M_bulk[k]}")
                    all_met = False
                    epsilon = epsilon_low[k] if np.abs(tolerances[k]) < 0.01 * M_bulk[k] else epsilon_high[k]
                    if tolerances[k] > 0:
                        partial_pressures.at[i, mol] += (-1) * epsilon
                    else:
                        partial_pressures.at[i, mol] += epsilon

            elif mol == 'Xe':       # For noble gases 

                if np.abs(tolerances[k]) > 1e-3 * M_bulk[k]:
                    print(f"Checking {mol}: Tolerance = {tolerances[k]}, Condition = {tolerances[k] > 0.1 * M_bulk[k]}")
                    all_met = False
                    epsilon = epsilon_low[k] if np.abs(tolerances[k]) < 0.01 * M_bulk[k] else epsilon_high[k]
                    if tolerances[k] > 0:
                        partial_pressures.at[i, mol] += (-1) * epsilon
                    else:
                        partial_pressures.at[i, mol] += epsilon

            # Break the loop if all tolerances are met
        if all_met:
            break

            # Liquid phase compositon : 

    for (mol,j) in zip(['H2O', 'CO2', 'NH3', 'CH4', 'Ar', 'Kr', 'Xe'], range(7)): 
        if mol == 'H2O':
            composition_liquidphase[mol][i] = res[j]
        else:
            composition_liquidphase[mol][i] = (res[j]) / (M_H2O*res[0])

    # composition_liquidphase['H2O'][i], composition_liquidphase['CO2'][i], composition_liquidphase['NH3'][i], composition_liquidphase['CH4'][i],\
    #     composition_liquidphase['Ar'][i], composition_liquidphase['Kr'][i], composition_liquidphase['Xe'][i] = res

    P = partial_pressures['H2O'][i] + partial_pressures['CO2'][i] + partial_pressures['NH3'][i] + partial_pressures['CH4'][i] \
        + partial_pressures['Ar'][i] + partial_pressures['Kr'][i] + partial_pressures['Xe'][i]

    y_H2O, y_CO2, y_NH3, y_CH4, y_Ar, y_Kr, y_Xe = partial_pressures['H2O'][i] / P, partial_pressures['CO2'][i] / P, partial_pressures['NH3'][i] / P, partial_pressures['CH4'] / P, \
        partial_pressures['Ar'][i] / P, partial_pressures['Kr'][i] / P, partial_pressures['Xe'][i] / P

    x_H2O, x_CO2, x_NH3, x_CH4, x_Ar, x_Kr, x_Xe = res

    P_CO2_beforechemical_eq = partial_pressures['CO2'][i]
    P_NH3_beforechemical_eq = partial_pressures['NH3'][i]

    #Compute chemial equlibrium to update PCO2 & PNH3

    species   = species[:9]
   # User input : total molality of CO2 and NH3 :

    m_CO2_tot = x_CO2 / (M_H2O*x_H2O)
    m_NH3_tot = x_NH3 / (M_H2O*x_H2O)

    # Initial guess to take from UNIQUAC

    # TKEL = ctypes.c_double(T)
    # FEED = (ctypes.c_double * 23)(55.5, m_NH3_tot, m_CO2_tot, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    #         # Call the function
    # BUBBLEP(ctypes.byref(TKEL), FEED, ctypes.byref(PBUB), PRES, ctypes.byref(NGAS),
    #             NBR, HGABS, ctypes.byref(AW), ctypes.byref(GM), ctypes.byref(ERROR))

    # # #          # Access the output values
    # # partial_pressures['H2O'][i] = PRES[0]*1e5
    # # partial_pressures['CO2'][i] = PRES[2]*1e5
    # # partial_pressures['NH3'][i] = PRES[1]*1e5

    # # # P0 = (PRES[0] + PRES[1] + PRES[2])*1e5    #[Pa]

    # m_CO2     = FEED[2]
    # m_NH3     = FEED[1]
    # m_HCO3m   = FEED[9]
    # m_Hp      = FEED[6]
    # m_NH4p    = FEED[5]
    # m_CO32m   = FEED[8]
    # m_OHm     = FEED[7]
    # m_NH2COOm = FEED[10]

    #z0 = np.array([0.99, m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm])
    # z0 = [0.99, 3e-4, 2.5e-1, 2.8e-1, 1.3e-9, 5.7e-1, 7.5e-2, 1e-5, 1.4e-1]

    # sol = least_squares(Chemical_eq_fun, z0, method='lm', ftol=1e-8, xtol=1e-8, verbose=1)

    # x_H2O, m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm = sol.x

    # x_HCO3m = m_HCO3m * (M_H2O*x_H2O)
    # x_Hp = m_Hp * (M_H2O*x_H2O)
    # x_NH4p = m_NH4p * (M_H2O*x_H2O)
    # x_CO32m = m_CO32m* (M_H2O*x_H2O)
    # x_OHm = m_OHm * (M_H2O*x_H2O)
    # x_NH2COOm = m_NH2COOm * (M_H2O*x_H2O)
    
    # z0 = [PsatH2O, 300, 30]

    # x_CO2 = m_CO2 * (M_H2O*x_H2O)
    # x_NH3 = m_NH3 * (M_H2O*x_H2O)

    # composition_liquidphase['H2O'][i], composition_liquidphase['CO2'][i], composition_liquidphase['NH3'][i] = x_H2O, x_CO2, x_NH3
    # composition_liquidphase['HCO3-'][i], composition_liquidphase['H+'][i], composition_liquidphase['NH4+'][i] = x_HCO3m, x_Hp, x_NH4p
    # composition_liquidphase['CO32-'][i], composition_liquidphase['OH-'][i], composition_liquidphase['NH2COO-'][i] = x_CO32m, x_OHm, x_NH2COOm

    # sol = least_squares(L_V_eq_fun, z0, method='lm', ftol=1e-8, xtol=1e-8, verbose=1)

    # partial_pressures['H2O'][i], partial_pressures['CO2'][i], partial_pressures['NH3'][i] = sol.x


    return composition_liquidphase.iloc[i], partial_pressures.iloc[i], P_CO2_beforechemical_eq, P_NH3_beforechemical_eq, M_atm

# %%
