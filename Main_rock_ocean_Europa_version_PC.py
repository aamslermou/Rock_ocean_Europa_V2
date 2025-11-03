## Main_rock_ocean_Europa.py

# %%                                     ##### Importing Libraries #####
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from scipy.optimize import fsolve

from PHREEQC_functions import create_pqi_file, run_phreeqc
database_file = "/Users/aamsler/Desktop/phreeqc-3.5.0-14000/database/core10.dat"  # Adjust the database path

from L_V_equilibrium_version_PC import L_V_equilibrium
from Saturation_pressure import P_sat_water
from Compo_ocean_Moon import compo_67P
# %%                                 ##### Constants and Parameters #####

G       = 6.67430e-11 #m3 kg−1 s−2
M_H2O   = 18.01528*1e-3 # [kg/mol]
rho_H2O = 1000      # [kg/m^3]
#rho_rock = 

N = 6 # Number of ocean boxes

# species is a panda dataframe containing the physical and chemical properties of the molecules of interest :
# formula : molecule's formula 
# Tc : critical temperature [K]
# Pc : critical pressure [Pa]
# w  : acentric factor []
# M  : molar mass [g/mol]
# charge : charge value 
# r_coeff : r parameter for gamma computation (Darde et al 2020)
# q_coeff : q parameter for gamma computation (Darde et al 2020) 
                                  
species = pd.DataFrame(np.array([['H2O', 647.3, 220.6*1e5, 0.3434, 18.01528, 0, 0.92, 1.40],
                         ['CO2', 304.19, 73.83*1e5, 0.224,  44.0095, 0, 0.75, 2.45],
                         ['NH3', 405.6, 112.8*1e5, 0.25, 17.03052, 0, 1.6292, 2.9852],
                         ['HCO3-', None, None, None, 61.0168, -1, 8.0756, 8.6806],
                         ['H+', None, None, None, 1.008, 1, 0.1378, 0.1e-15],
                         ['NH4+', None, None, None, 18.03846, 1, 4.8154, 4.6028],
                         ['CO32-', None, None, None, 60.0089, -2, 10.828, 10.769],
                         ['OH-', None, None, None, 17.00734, -1, 9.3973, 8.8171],
                         ['NH2COO-', None, None, None, 108.9936, -1, 4.3022, 4.1348],
                         ['CH4', 190.4, 46.0*1e5, 0.011, 16.0425, 0, None, None],
                         ['Ar', 150.8, 48.7*1e5, 0.001, 39.948, 0, None, None],
                         ['Kr', 209.4, 55.0*1e5, 0.005, 83.798, 0, None, None],
                         ['Xe', 289.7, 58.4*1e5, 0.008, 131.29, 0, None, None]]),
                  
              columns = ['formula', 'Tc', 'Pc', 'w', 'M','charge', 'r_coeff', 
                'q_coeff'])

species_phi   = pd.concat([species[:3],species[9:]], ignore_index=True)
                                        #######


# k_coeff is a panda dataframe containing the binary interaction parameters kij, aiming to correct A value in Peng-Robinson equation of state
# Such parameters are symetrical, and ajusted on experimental data 
# source : cf Amsler Moulanier et al. (2025a),                                                                     
k_coeff = np.array([[0,0.1896,0,0,0,0,0],       # H2O
                     [0.1896,0,0,0.084,0,0,0],  # CO2
                     [0,0,0,0,0,0,0],           # NH3
                     [0,0.084,0,0,0,0,0],       # CH4
                     [0,0,0,0,0,0,0],           # Ar
                     [0,0,0,0,0,0,0],           # Kr
                     [0,0,0,0,0,0,0]])          # Xe

                                        ########
u0_coeff = pd.DataFrame(np.array([[0, 8.8383,594.72,577.05,10000,52.7305,361.39,600.50,28.2779],
                                 [8.8383,302.25,2500.0,526.305,1e9,-424.01,2500.0,2500.0,2500.0],
                                 [594.72,2500,1090.8,534.01,1e9,785.98,524.13,1733.9,498.15],
                                 [577.05,526.305,534.01,771.04,1e9,505.55,800.01,2500.0,613.25],
                                 [10000,1e9,1e9,1e9,0,1e9,1e9,1e9,1e9],
                                 [52.7305,-424.01,785.98,505.55,1e9,0,226.60,1877.9,44.849],
                                 [361.39,2500.0,524.13,800.01,1e9,226.6,1458.3,1588,2500.0],
                                 [600.5,2500.0,1733.9,2500,1e9,1877.9,1588,1562.9,2500.0],
                                 [28.2779, 2500.0, 498.15, 613.25, 1e9, 44.849, 2500.0, 2500.0, 3343.1]]),
                       
                        columns = ['H2O', 'CO2', 'NH3', 'HCO3-', 'H+','NH4+','CO32-','OH-','NH2COO-'])

                                        ########
                                    
# ut_coeff is a panda dataframe containing the ut parameters involved in the activity coefficient compuation
# source : Darde et al 2010
ut_coeff = pd.DataFrame(np.array([[0,0.86293,7.1827,-0.38795,0,0.50922,3.3516,8.5455,8.0238],
                                 [0.8629,0.35870,0,-3.7340,0,8.6951,0,0,0],
                                 [7.1827,0,7.0912,5.3111,0,6.1271,4.9305,0.1364,6.6532],
                                 [-0.388,-3.734,5.3111,-0.01981,0,-0.00795,1.7241,0,3.0580],
                                 [0,0,0,0,0,0,0,0,0],
                                 [0.5092,8.6951,6.1227,-0.00795,0,0,4.0555,0.34921,12.047],
                                 [3.3516,0,4.9305,1.724,0,4.056,-1.3448,2.7496,0],
                                 [8.5455,0,0.1364,0,0,0.3492,2.7496,5.6169,0],
                                 [8.0238,0,6.6532,3.0580,0,12.047,0,0,-15.920]]),
                       
                        columns = ['H2O', 'CO2', 'NH3','HCO3-','H+','NH4+','CO32-','OH-','NH2COO-'])

# Molar masses in kg/mol
M  = pd.DataFrame(np.array([[18.01528e-3, 44.0095e-3, 16.0425e-3, 39.948e-3, 83.798e-3, 131.29e-3, 17.03052e-3, 31.9988e-3, 28.0101E-3, 28.0134e-3, 61.0168e-3, 18.03846e-3, 60.0089e-3, 108.9936e-3]]),      # Molar masses in [kg/mol]

                                    columns= ['H2O', 'CO2', 'CH4', 'Ar', 'Kr', 'Xe', 'NH3', 'O2', 'CO', 'N2', 'HCO3-', 'NH4+', 'CO32-', 'NH2COO-'])

#name_species = np.array(['H2O', 'CO2', 'NH3', 'CH4', 'CO', 'O2', 'N2', 'Ar', 'Kr', 'Xe'])
name_species = np.array(['H2O', 'CO2', 'NH3', 'CH4', 'Ar', 'Kr', 'Xe'])
molar_masses = np.array([float(M[mol]) for mol in name_species])
compo = np.array([float(compo_67P[mol]) for mol in name_species])


# PHRREQ-C
database_file = "/Users/aamsler/Desktop/phreeqc-3.5.0-14000/database/core10.dat"  # Adjust the database path

# %%                                    ##### Initialization of the model #####

# Yannis' input : 

Accretion_variables = pd.read_csv('1Myr_270K.csv', delimiter = ";")

time        = Accretion_variables['t']
M_evolution = Accretion_variables['M_Europa']
T_evolution = Accretion_variables['T_surf']
R_evolution = Accretion_variables['R_Europa']

nsteps      = int(len(time)/10)

M_impactors = M_evolution[9] - M_evolution[0]       # Mass of impactors contributing to accretion process
delta_time  = time[1] - time[0]

M_bulk = pd.DataFrame(np.zeros((len(time),7)), 
                      columns=['H2O', 'CO2', 'NH3', 'CH4', 'Ar', 'Kr', 'Xe'])

rho_rock = 3000 # [kg/m^3]
r_rock   = 1    # [m]                    [between 1 cm and 100m]
eta      = 0.0017914 # [Pa.s]  # Viscosity of the ocean  
# %%                             ##### Definition of the physical system #####

## Rocky core, bottom box
rocky_core = pd.DataFrame(np.array([np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps)]).T, columns=['T', 'P', 'W:R', 'V_H2O', 'V_rock'])
output_dir = "Rocky_core"
os.makedirs(output_dir, exist_ok=True)
## Ocean, N boxes
box_data = {}

for n in range(N):
    box_data[f"B_{n}_{n+1}"] = pd.DataFrame(np.array([np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps)]).T, columns=['T', 'P', 'Depth', 'W:R', 'CO2_tot', 'NH3_tot'])
    output_dir = f"Box_{n}_{n+1}"
    os.makedirs(output_dir, exist_ok=True)

composition_liquidphase = pd.DataFrame(np.zeros((nsteps,14)),

                                    columns= ['H2O', 'CO2', 'NH3', 'HCO3-', 'H+', 'NH4+', 'CO32-', 'OH-', 'NH2COO-', 'CH4', 'Ar', 'Kr', 'Xe', 'pH'])

partial_pressures       =  pd.DataFrame(np.zeros((nsteps,7)),

                                    columns= ['H2O', 'CO2', 'NH3', 'CH4', 'Ar', 'Kr', 'Xe'])

# %%                                    ##### Definition of the functions #####
 
# Stokes flow
def stokes_flow(eta, d, g, r_rock, rho_rock, rho_water):
    return ( (9/2) * (eta*d) / (g * (rho_rock - rho_water) * r_rock**2) )

# Hydrostatic pressure
def hydrostatic_pressure(depth, density, g):
    return density * g * depth

# %%
# Computation of the atmosphere and upper-ocean composition at t = 0

#T          = T_evolution[0]                      # Surface temperature at the start of accretion
T = 273.15      #We start the simulation at 0°C

M_Europa   = M_evolution[0]            # [kg]

M_bulk_tot = 0.08 * M_Europa  # 8% water content

R_Europa   = R_evolution[0]               # [m]
g_Europa   = G * M_Europa / (R_Europa**2) # [m/s^2]    

# Computation of the ocean's depth at t = 0
M_tot    = np.sum(molar_masses*compo)
M_bulk.iloc[0] = compo * M_bulk_tot * ( molar_masses / M_tot )
M_ocean  = M_bulk['H2O'][0] - ( (4 * np.pi * R_Europa**2 * P_sat_water(T)) / g_Europa)
V_Ocean  = M_ocean / rho_H2O
d_Ocean  = R_Europa - (R_Europa**3 - 3 *V_Ocean / (4*np.pi)) ** (1/3)

# Computation of the rocky core layer thickness at t = 0
R_core   = R_Europa - d_Ocean
d_porous = 0.1 * R_core                                        # Porous layer thickness
V_core   = 4/3 * np.pi * (R_core**3 - (R_core - d_porous)**3)   # Volume of the porous layer
W_R_core = 1                                                  # Initial W/R ratio of the rocky core layer

# Volume of H2O pushed out of the porous layer due to compaction : 
V_H2O_out = (0.45 * M_impactors ) / rho_rock    # Volume of water out is equal to the volume of rock falling at the bottom of the sea. 
# Compaction percentage : 
Compac = np.linspace(0.1, 0.6, nsteps)  # Compaction percentage increasing linearly from 10% to 40% over the simulation time, as pressure at the bottom increases

f1 = lambda W_R_core, V_H2O, V_rock, rho_rock, rho_H2O: W_R_core - (V_H2O / V_rock) * (rho_H2O / rho_rock)
f2 = lambda V_H2O, V_rock, R_Europa, d_Ocean, d_porous: V_H2O - ( (4/3 * np.pi * (R_Europa - d_Ocean)**3 - 4/3 * np.pi * (R_Europa - d_Ocean - d_porous)**3) - V_rock)

def fun_to_solve(z):

    V_H2O, V_rock = z

    return [f1(W_R_core, V_H2O, V_rock, rho_rock, rho_H2O), f2(V_H2O, V_rock, R_Europa, d_Ocean, d_porous)]

rocky_core['V_H2O'][0], rocky_core['V_rock'][0] = fsolve(fun_to_solve, np.array([5e13, 5e13]))

W_R_core_range = np.linspace(1, 0.1, nsteps)  # W/R ratio of the rocky core decreasing linearly from 1 to 0.1 over the simulation time
# In the case t_stockes is too fast !!!

    ## Computation of the atmosphere-ocean interface equilibrium :

#initial_guess = np.array([29.6e5, 300, 4e5, 1782, 304, 224])
initial_guess = np.array([3.22e4, 140, 844, 2.5, 0.5, 0.5]) # Initial guess for the partial pressures 
j = 0
composition_liquidphase.iloc[j], partial_pressures.iloc[j], PCO2_beforechem, PNH3_beforechem, M_atm = L_V_equilibrium(T, name_species, M_bulk.iloc[j], initial_guess, composition_liquidphase, partial_pressures, j, d_Ocean, R_Europa, g_Europa, u0_coeff, ut_coeff, k_coeff, species_phi, species)

    ## Computation of geochemical equilibrium in each box at t = 0

# %%
# for n in range(0,N):
    
#     output_dir = './' + f"Box_{n}_{n+1}" +'/'
#     depth      = d_Ocean * (n + 1/2) / (N)        # Mean depth of the box
#     V_box      = 4/3 * np.pi * ((R_Europa - d_Ocean * n/N)**3 - (R_Europa - d_Ocean * (n+1)/N)**3)
#     Water_mass_box = V_box * rho_H2O
#     R = Water_mass_box/10    # Rock mass in ocean at t=0
#     box_data[f"B_{n}_{n+1}"]['P'][i]   = hydrostatic_pressure(depth, rho_H2O, g_Europa)*1e-5 # Hydrostatic pressure in bar
#     #box_data[f"B_{n}_{n+1}"]['P'][i] = 50 # Hydrostatic pressure in Pa
#     box_data[f"B_{n}_{n+1}"]['W:R'][i] = Water_mass_box/R - n
#     if T==273.15:
#         box_data[f"B_{n}_{n+1}"]['T'][i] = 0.01
#     else:
#         box_data[f"B_{n}_{n+1}"]['T'][i]   = int(T - 273.15) # Temperature in Celsius
#     # Change the input files and compute the new geochemical equilibrium with PHREEQC

#     input_file = "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".pqi"
#     data_file  = "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".dat"
#     create_pqi_file("Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".pqi", box_data[f"B_{n}_{n+1}"]['W:R'][i], box_data[f"B_{n}_{n+1}"]['T'][i], box_data[f"B_{n}_{n+1}"]['P'][i],output_dir,data_file)
#     outputfile = "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".txt"
#     run_phreeqc(input_file, outputfile, output_dir, database_file)
#     new_input_path = os.path.join(output_dir, input_file)
#     os.rename(input_file, new_input_path)  # Move the input file to the output directory

# rocky_core['P'][i] = hydrostatic_pressure(d_Ocean, rho_rock, g_Europa)*1e-5 # Hydrostatic pressure
# if T==273.15:
#     rocky_core['T'][i] = 0.01
# else:
#     rocky_core['T'][i] = T  # Update temperature

# rocky_core['W:R'][i] = W_R_core # Initial W/R ratio at the porous interface between rocky mantle and ocean
# output_dir = './' + "Rocky_core" +'/'
# input_file = "Rocky_core" + "_t_" + str(int(time[i])) + ".pqi"
# data_file  = "Rocky_core" + "_t_" + str(int(time[i])) + ".dat"
# create_pqi_file( "Rocky_core" + "_t_" + str(int(time[i])) + ".pqi",rocky_core['W:R'][i], rocky_core['T'][i], rocky_core['P'][i], output_dir,data_file)
# outputfile = "Rocky_core" + "_t_" + str(int(time[i])) + ".txt"
# run_phreeqc(input_file, outputfile, output_dir, database_file)
# new_input_path = os.path.join(output_dir, input_file)
# os.rename(input_file, new_input_path)  # Move the input file to the output directory
# %%                                   ##### Reset box, values at t0 of Mbulk and Matm #####

# M_bulk = np.array([7.68925202e+17, 8.82850915e+16, 4.87018872e+15, 2.32805933e+15,
#        9.88931279e+12, 1.75255967e+12, 1.34488754e+12])
# M_atm = np.array([8.86787339e+16, 3.90544449e+14, 2.32437427e+15, 7.56247840e+12,
#        1.74328070e+12, 1.34395100e+12])
# %%                                      ##### Main loop of the model #####

d_Ocean_list = []

for i in range(1, len(time), 10):
    t = time[i]
    j += 1

    M_ice_new = M_impactors * 0.08   # 8% of the ice contained in the impactors contributes to the ocean formation   
    #M_ice_new = M_impactors * 0.55  # Impactors are made of 55% of ice (Callisto ice fraction)
    M_atm    += M_ice_new
    M_bulk.iloc[j] = M_bulk.iloc[j-1] + compo * M_ice_new * (molar_masses / M_tot)

    factor = M_bulk['H2O'][j] / M_bulk['H2O'][j]  # Increase factor for the bulk composition

    T_surface = T_evolution[i]
        
    # Computation of Europa's characteristics and its the ocean's depth at time t

    M_Europa = M_evolution[i]            # [kg]
    R_Europa = R_evolution[i]            # [m]
    g_Europa = G * (M_Europa / R_Europa**2)     

    P_satH2O = P_sat_water(T_surface)
    M_H2O_atm_max = (4 * np.pi * R_Europa**2 * P_satH2O) / g_Europa

    M_ocean = M_bulk['H2O'][j] - M_H2O_atm_max
    V_Ocean = M_ocean / rho_H2O

    d_Ocean = R_Europa - (R_Europa**3 - 3 * V_Ocean / (4*np.pi)) ** (1/3)
    d_Ocean_list.append(d_Ocean)
    
    # Computation of the rocky core layer thickness at time t

    # R_core   = R_Europa - d_Ocean
    # d_porous = 0.01 * R_core                                        # Porous layer thickness (1% of the mantle)
    # V_core   = 4/3 * np.pi * (R_core**3 - (R_core - d_porous)**3)   # Volume of the porous layer

    # Water_mass_box = d_Ocean / N          # Water mass in each box

    # L-V equilibrium at interface     

    # if i<20:
    #     initial_guess = np.array([PCO2_beforechem * (factor * 0.94),
    #                             PNH3_beforechem * (factor * 0.85),
    #                             partial_pressures['CH4'][i-1] * (factor*0.9325) ,
    #                             partial_pressures['Ar'][i-1] * factor,
    #                             partial_pressures['Kr'][i-1] * factor,
    #                             partial_pressures['Xe'][i-1] * factor])
    # else:
    initial_guess = np.array([PCO2_beforechem * (factor*0.95),
                                PNH3_beforechem * (factor * 0.92),
                                partial_pressures['CH4'][j-1] * (factor*0.98) ,
                                partial_pressures['Ar'][j-1] * factor,
                                partial_pressures['Kr'][j-1] * factor,
                                partial_pressures['Xe'][j-1] * factor])
        
    # initial_guess = np.array([PCO2_beforechem * (factor[1]),
    #                           PNH3_beforechem * (factor[2]),
    #                           partial_pressures['CH4'][i-1] * (factor[3]) ,
    #                           partial_pressures['Ar'][i-1] * factor[4],
    #                           partial_pressures['Kr'][i-1] * factor[5],
    #                           partial_pressures['Xe'][i-1] * factor[6]])      
     
    composition_liquidphase.iloc[j], partial_pressures.iloc[j], PCO2_beforechem, PNH3_beforechem, M_atm = L_V_equilibrium(T_surface, name_species, M_bulk.iloc[j], initial_guess, composition_liquidphase, partial_pressures, j, d_Ocean, R_Europa, g_Europa, u0_coeff, ut_coeff, k_coeff, species_phi, species)

    print(i/len(time)*100,'%')

    composition_liquidphase.to_csv('composition_liquidphase_1Myr_270.csv', index=False)
    partial_pressures.to_csv('partial_pressures_1Myr_270K.csv', index=False)
    #Stokes_time_scale = stokes_flow(eta,d_Ocean, g_Europa, r_rock, rho_rock, rho_H2O)
    #if Stokes_time_scale > delta_time:
        
    # for n in range(0,N):

    #     output_dir = './' + f"Box_{n}_{n+1}" +'/'
    #     depth = d_Ocean * (n + 1/2) / (N)        # Mean depth of the box
    #     box_data[f"B_{n}_{n+1}"]['P'][i] = hydrostatic_pressure(depth, rho_rock, g_Europa)*1e-5 # Hydrostatic pressure
    #     #dR_bottom = dR_up
    #     #dR_up = function for sedimentation at depth n (dR >0)
    #     #R = eval(f"B_{n}_{n+1}")['W:R'][i-1] + dR_up - dR_bottom
    #     #eval(f"B_{n}_{n+1}")['W:R'][i] = Water_mass_box/R
    #     box_data[f"B_{n}_{n+1}"]['W:R'][i] = 10
    #     # Change the input files and compute the new geochemical equilibrium with PHREEQC

    #     if T==273.15:
    #         box_data[f"B_{n}_{n+1}"]['T'][i] = 0.01
    #     else:
    #         box_data[f"B_{n}_{n+1}"]['T'][i] = int(T_surface - 273.15) # Temperature in Celsius

    #     input_file = "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".pqi"
    #     data_file  = "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".dat"
    #     create_pqi_file("Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".pqi", box_data[f"B_{n}_{n+1}"]['W:R'][i], box_data[f"B_{n}_{n+1}"]['T'][i], box_data[f"B_{n}_{n+1}"]['P'][i],output_dir,data_file)
    #     outputfile = "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(int(time[i])) + ".txt"
    #     run_phreeqc(input_file, outputfile, output_dir, database_file)
    #     new_input_path = os.path.join(output_dir, input_file)
    #     os.rename(input_file, new_input_path)  # Move the input file to the output directory
       
    #     # dW_core = #function for compaction (dW >0)
    #     # #dW_top = M_ice_new
    #     # W = rocky_core['W:R'][i]*Rock_mass_core - dW
    #     # rocky_core['W:R'][i] = W/Rock_mass_core     # Update W/R, it should decrease because of compaction
    #     # rocky_core['T'][i] = rocky_core['T'][t-1] + dT #Update temperature

    #     # Compaction of the porous rocky layer : computation of the new W/R ratio
    
    # #W_R_core = rocky_core['W:R'][i-1]  # Previous W/R ratio
    # W_R_core = 1
    # V_H2O, V_rock = fsolve(fun_to_solve, np.array([rocky_core['V_H2O'][i-1], rocky_core['V_rock'][i-1]]))

    # V_H2O_new = V_H2O - V_H2O_out * Compac[i] # We assume that a certain % of the water is pushed out of the porous layer 

    # if V_H2O_new < 0:
    #     print('error: no more water in the porous layer')
    #     break

    # rocky_core['V_H2O'][i] = V_H2O_new
    # rocky_core['V_rock'][i] = V_rock

    # rocky_core['W:R'][i] = (V_H2O_new / V_rock) * (rho_H2O / rho_rock)  # Update W/R, it should decrease because of compaction
    # # rocky_core['W:R'][i] = W_R_core_range[i]

    # rocky_core['P'][i] = hydrostatic_pressure(d_Ocean, rho_rock, g_Europa)*1e-5 # Hydrostatic pressure
    # if T==273.15:
    #     rocky_core['T'][i] = 0.01
    # else:
    #     rocky_core['T'][i] = T_surface  # Update temperature

    # output_dir = './' + "Rocky_core" +'/'
    # input_file = "Rocky_core" + "_t_" + str(int(time[i])) + ".pqi"
    # data_file  = "Rocky_core" + "_t_" + str(int(time[i])) + ".dat"
    # create_pqi_file( "Rocky_core" + "_t_" + str(int(time[i])) + ".pqi",rocky_core['W:R'][i], rocky_core['T'][i], rocky_core['P'][i], output_dir,data_file)
    # outputfile = "Rocky_core" + "_t_" + str(int(time[i])) + ".txt"
    # run_phreeqc(input_file, outputfile, output_dir, database_file)
    # new_input_path = os.path.join(output_dir, input_file)
    # os.rename(input_file, new_input_path)  # Move the input file to the output directory
# %%
# composition_liquidphase.to_csv('composition_liquidphase_1Myr_270.csv', index=False)
# partial_pressures.to_csv('partial_pressures_1Myr_270K.csv', index=False)
# %%
