# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sszpalette

from labellines import labelLine, labelLines
from matplotlib import gridspec

# %%
N = 6

pH = np.zeros(N) 
m_CO2  = np.zeros(N)
m_NH3  = np.zeros(N)
m_CH4  = np.zeros(N)
K = np.zeros(N)
Cl = np.zeros(N)
C = np.zeros(N)
Na = np.zeros(N)

d_ocean = 125 # [km]
depth = np.array([d_ocean * (n + 1/2) / N for n in range(N)]) # [km]

colorsmaps = sszpalette.register()
colorsmaps

cmap = 'contrasting12'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), layout='constrained')

for n in range(N):
    directory = './' + f"Box_{n}_{n+1}" +'/'
    data = pd.read_csv(directory + "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(0) + '.dat', delimiter = "\t")
    data = pd.DataFrame(data)
    pH[n]    = data['          pH'][1]
    K[n]     = data['           K'][1]
    Cl[n]    = data['          Cl'][1]
    C[n]     = data['           C'][1]
    Na[n]    = data['          Na'][1]
    #m_CO2[n] = data['       m_CO2'][1]
#m_NH3[n] = data[       'm_NH3'][1]
    #m_CH4[n] = data['       m_CH4'][1]


plt.set_cmap('contrasting12')
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.get_cmap().colors)
ax1.plot(pH, depth, label = 'pH', color = 'C0')
#ax2.plot(m_CO2, depth, label = 'm_CO2')
#ax1.plot(m_NH3, depth, label = 'm_NH3')

ax1.invert_yaxis()  # Invert the x-axis to show decreasing temperature

#ax2.plot(m_CH4, depth, label = 'm_CH4')
ax2.plot(K, depth, label = 'K', color = 'C1')
ax2.plot(Cl, depth, label = 'Cl', color = 'C2')
ax2.plot(C, depth, label = 'C', color = 'C3')
ax2.plot(Na, depth, label = 'Na', color = 'C4')

labelLines(ax2.get_lines(), zorder=2.5)

ax1.set_xlabel('pH')
ax2.set_xlabel('Concentration [mol/kg]')
ax1.set_ylabel('Depth [km]')

ax2.set_ylabel('Depth [km]')
ax2.invert_yaxis()  # Invert the x-axis to show decreasing temperature
ax2.set_xscale('log')

plt.show()
# %% Partial pressure vs time

# Plot partial pressures from CSV file
pp_data = pd.read_csv('Output/partial_pressures_1Myr_270K.csv')
comp_data = pd.read_csv('Output/composition_liquidphase_1Myr_270.csv')

time_plot = np.linspace(100*1e-3, 1e6*1e-3, 1000)
molecule_names = ['H2O', 'CO2', 'NH3','CH4']

fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12, 6), layout='constrained')
for molecule in molecule_names:
    if molecule=='H2O':
        ax1.plot(time_plot, pp_data[molecule], label=molecule, color = 'lightblue')
    else:
        ax1.plot(time_plot, pp_data[molecule], label=molecule)

labelLines(ax1.get_lines(), zorder=2.5, xvals=[180,300,190,400])

molecule_names = ['CO2', 'NH3', 'HCO3-', 'NH4+', 'CH4']

for molecule in molecule_names:
    if molecule=='HCO3-':
        ax2.plot(time_plot, comp_data[molecule], label=molecule, color = 'cyan')
    elif molecule=='NH4+':
        ax2.plot(time_plot, comp_data[molecule], label=molecule, color = 'black')
    else:
        ax2.plot(time_plot, comp_data[molecule], label=molecule)

labelLines(ax2.get_lines(), zorder=2.5, xvals=[800,160,260,454,300])

ax1.set_xlabel('Time [kyr]', fontsize=14)
ax1.set_ylabel('Partial Pressure [Pa]', fontsize=14)
ax1.set_yscale('log')
ax1.set_title('Primordial atmosphere', fontsize=15)

ax2.set_xlabel('Time [kyr]', fontsize=14)
ax2.set_ylabel('Molality [mol/kg]', fontsize=14)
ax2.set_yscale('log')
ax2.set_title('Primordial ocean', fontsize=15)
plt.savefig('Europa_evol_1Myr_270K.svg', format='svg')
plt.show()

# %% Partial pressure vs time zoom in

# Plot partial pressures from CSV file
pp_data = pd.read_csv('partial_pressures_1Myr_270K.csv')
pp_data_zoom = pp_data.iloc[0:25]

time_plot = np.linspace(100, 25*1000, 25)
molecule_names = ['H2O', 'CO2', 'NH3','CH4']
fig, ax1 = plt.subplots(layout='constrained')
for molecule in molecule_names:
    if molecule=='H2O':
        ax1.plot(time_plot, pp_data_zoom[molecule], label=molecule, color = 'lightblue')
    else:
        ax1.plot(time_plot, pp_data_zoom[molecule], label=molecule)
labelLines(ax1.get_lines(), zorder=2.5, xvals=[2,3,2.5,4])
plt.xlabel('Time [kyr]', fontsize=14)
plt.ylabel('Partial Pressure [Pa]', fontsize=14)
plt.yscale('log')
#plt.savefig('Europa_evol_atm_10kyr_270K.svg', format='svg')
plt.show()
# %% LV interface composition vs time

# Plot composition of the liquid phase from CSV file
comp_data = pd.read_csv('composition_liquidphase_1Myr_270.csv')

time_plot = np.linspace(100*1e-3, 1e6*1e-3, 1000)
molecule_names = ['CO2', 'NH3', 'HCO3-', 'NH4+', 'CH4']

fig, ax1 = plt.subplots(layout='constrained')
for molecule in molecule_names:
    if molecule=='HCO3-':
        ax1.plot(time_plot, comp_data[molecule], label=molecule, color = 'cyan')
    elif molecule=='NH4+':
        ax1.plot(time_plot, comp_data[molecule], label=molecule, color = 'black')
    else:
        ax1.plot(time_plot, comp_data[molecule], label=molecule)

labelLines(ax1.get_lines(), zorder=2.5, xvals=[800,160,260,454,300])

plt.xlabel('Time [kyr]', fontsize=14)
plt.ylabel('Molality [mol/kg]', fontsize=14)
plt.yscale('log')
plt.savefig('Europa_evol_liquid_1Myr_270K.svg', format='svg')
plt.show()
# %% Rocky core vs time

time_plot = np.linspace(200*1e-3,998200*1e-3,998)
time_file = range(200,998200,1000)

N = len(time_file)
pH = np.zeros(N)
m_CO2  = np.zeros(N)
m_NH3  = np.zeros(N)
m_CH4  = np.zeros(N)
K = np.zeros(N)
Cl = np.zeros(N)
C = np.zeros(N)
Na = np.zeros(N)
Ca = np.zeros(N)
Mg = np.zeros(N)
Fe = np.zeros(N)
Si = np.zeros(N)
HCO3m = np.zeros(N)
NaCl = np.zeros(N)
C_m4 = np.zeros(N)
Cronstedtite = np.zeros(N)
Lizardite = np.zeros(N)
Magnetite = np.zeros(N)
Hydroxyapatite = np.zeros(N)
Pyrite = np.zeros(N)    
Magnesite = np.zeros(N)
Hematite = np.zeros(N)
Dolomite = np.zeros(N)
Calcite = np.zeros(N)
Talc = np.zeros(N)
Chamosite = np.zeros(N)
Saponite = np.zeros(N)

for n in range(N):
    directory = 'Output/' + "Rocky_core" +'/'
    data = pd.read_csv(directory + "Rocky_core" + "_t_" + str(time_file[n]) + '.dat', delimiter = "\t")
    data = pd.DataFrame(data)
    pH[n]    = data['          pH'][1]
    K[n]     = data['           K'][1]
    Cl[n]    = data['          Cl'][1]
    C[n]     = data['           C'][1]
    Na[n]    = data['          Na'][1]
    Ca[n]    = data['          Ca'][1]
    Mg[n]    = data['          Mg'][1]
    Fe[n]    = data['          Fe'][1]
    Si[n]    = data['          Si'][1]
    HCO3m[n] = data['     m_HCO3-'][1]
    m_NH3[n] = data['       m_NH3'][1]
    m_CO2[n] = data['       m_CO2'][1]
    C_m4[n]  = data['        C(4)'][1]
    NaCl[n]         =  data['      m_NaCl'][1]
    Cronstedtite[n] = data['Cronstedtite-7A'][1]
    Lizardite[n]  = data['   Lizardite'][1]
    Magnetite[n]  = data['   Magnetite'][1]
    Hydroxyapatite[n] = data['Hydroxyapatite'][1]
    Pyrite[n]     = data['      Pyrite'][1]
    Magnesite[n]  = data['   Magnesite'][1]
    Hematite[n]   = data['    Hematite'][1]
    Dolomite[n]   = data['    Dolomite'][1]
    Talc[n]       = data['        Talc'][1]
    Calcite[n]    = data['     Calcite'][1]
    Chamosite[n]  = data['   Chamosite'][1]
    Saponite[n]   = data['Saponite-Fe-Na'][1]

colorsmaps = sszpalette.register()
colorsmaps

cmap = 'contrasting12'

#from matplotlib import gridspec
# Use gridspec to set height ratios: ax1 is half the height of ax2 and ax3

# fig = plt.figure(figsize=(5, 11), layout='constrained')
# gs = gridspec.GridSpec(3, 1, height_ratios=[1, 2, 2], figure=fig)

# ax1 = fig.add_subplot(gs[0])
# ax2 = fig.add_subplot(gs[1])
# ax3 = fig.add_subplot(gs[2])
fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(12, 6), layout='constrained')
plt.set_cmap('contrasting12')
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.get_cmap().colors)
ax1.plot(time_plot, pH, label='pH', color='C0')

ax2.plot(time_plot, K, label='K', color='C1')
ax2.plot(time_plot, Cl, label='Cl', color='C2')
ax2.plot(time_plot, C, label='C', color='C3')
ax2.plot(time_plot, C_m4, label='C(4)', color='C19')

ax2.plot(time_plot, Na, label='Na', color='C4')
ax2.plot(time_plot, Ca, label='Ca', color='C5')
#ax2.plot(time_plot, HCO3m, label='HCO3-', color='C9')
ax2.plot(time_plot, Mg, label='Mg', color='C6')
ax2.plot(time_plot, Fe, label='Fe', color='C7')
ax2.plot(time_plot, Si, label='Si', color='C8')
#ax2.plot(time_plot, NaCl, label='NaCl', color='C10')

labelLines(ax2.get_lines(), zorder=2.5)

ax1.set_ylabel('pH', fontsize=14)
ax2.set_ylabel('Molality [mol/kg]', fontsize=14)
ax2.set_yscale('log')
ax2.set_xlabel('time [kyr]', fontsize=14)
ax1.set_xlabel('time [kyr]', fontsize=14)

ax3.plot(time_plot, Calcite, label='Calcite', color='C8')
ax3.plot(time_plot, Cronstedtite, label='Cronstedtite', color='C1')
ax3.plot(time_plot, Lizardite, label='Lizardite', color='C2')
ax3.plot(time_plot, Saponite, label='Saponite-Fe-Na', color='C3')
ax3.plot(time_plot, Chamosite, label='Chamosite', color='C9')
#ax3.plot(time_plot, Hydroxyapatite, label='Hydroxyapatite', color='C4')
ax3.plot(time_plot, Pyrite, label='Pyrite', color='C5')
ax3.set_yscale('log')
ax3.set_ylabel('Mineral amount [mol/kg]', fontsize=14)
ax3.set_xlabel('time [kyr]', fontsize=14)

labelLines(ax3.get_lines(), zorder=2.5, xvals=[50, 400, 500, 600, 200, 400])
#plt.savefig('Rocky_mantle_1Myr_270K.pdf', format='pdf')
plt.show()
# %% Composition of the water column at a given time 

# a t=0
N = 10

pH = np.zeros(N) 
pe = np.zeros(N)
m_CO2  = np.zeros(N)
m_NH3  = np.zeros(N)
m_CH4  = np.zeros(N)
K = np.zeros(N)
Cl = np.zeros(N)
C_4 = np.zeros(N)
C_m4 = np.zeros(N)
N_m3 = np.zeros(N)
Na = np.zeros(N)
Ca = np.zeros(N)
Mg = np.zeros(N)
Fe = np.zeros(N)
Calcite = np.zeros(N)
Cronstedtite = np.zeros(N)
Lizardite = np.zeros(N)
Saponite = np.zeros(N)
Chamosite = np.zeros(N)
Pyrite = np.zeros(N)

d_ocean = 5084
depth = np.array([d_ocean * (n + 1/2) / N for n in range(N)]) # [km]

for n in range(N):
    directory = 'Output/' + f"Box_{n}_{n+1}" +'/'
    data = pd.read_csv(directory + "Box_" + str(n) + "_" + str(n+1) + "_t_" + str(0) + '.csv', delimiter = ",")
    data = pd.DataFrame(data)
    pH[n]    = data['pH'][0]
    pe[n]    = data['pe'][0]
    K[n]     = data['m_K'][0]
    Cl[n]    = data['m_Cl'][0]
    C_4[n]     = data['m_C(4)'][0]
    C_m4[n]    = data['m_C(-4)'][0]
    N_m3[n]    = data['m_N(-3)'][0]
    Na[n]    = data['m_Na'][0]
    Ca[n]    = data['m_Ca'][0]
    Mg[n]    = data['m_Mg'][0]
    Fe[n]    = data['m_Fe'][0]
    Calcite[n] = data['calcite'][0]
    Cronstedtite[n] = data['cronstedtite'][0]
    Lizardite[n]  = data['lizardite'][0]
    Saponite[n]   = data['saponite-Fe-Na'][0]
    Chamosite[n]  = data['chamosite'][0]
    Pyrite[n]     = data['pyrite'][0]

colorsmaps = sszpalette.register()
colorsmaps

cmap = 'contrasting12'

fig, (ax1,ax2,ax3) = plt.subplots( 1,3, figsize=(12, 6), layout='constrained')

plt.set_cmap('contrasting12')
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.get_cmap().colors)
ax1.plot(pH, depth, label='pH', color='C0')

ax2.plot(K,  depth,label='K', color='C1')
ax2.plot(Cl, depth, label='Cl', color='C2')
ax2.plot(C_4,  depth,label='C(4)', color='C3')
ax2.plot(C_m4, depth, label='C(-4)', color='C19')
ax2.plot(N_m3, depth, label='N(-3)', color='C20')
ax2.plot(Na, depth, label='Na', color='C4')
ax2.plot(Ca, depth, label='Ca', color='C5')
ax2.plot(Mg, depth, label='Mg', color='C6')
ax2.plot(Fe, depth, label='Fe', color='C7')

labelLines(ax2.get_lines(), zorder=2.5, xvals=[1e-6,2e-3,3e-3,4e-3,5e-3,6e-6,7e-3,8e-3,9e-3])
#ax2.legend(fontsize=8, loc='upper left')
ax1.set_xlabel('pH', fontsize=14)
ax2.set_xlabel('Molality [mol/kg]', fontsize=14)
ax2.set_xscale('log')
ax1.set_ylabel('Depth [m]', fontsize=14)
ax2.set_ylabel('Depth [m]', fontsize=14)

ax3.plot(Calcite, depth, label='Calcite', color='C8')
ax3.plot(Cronstedtite, depth, label='Cronstedtite', color='C1')
ax3.plot(Lizardite, depth, label='Lizardite', color='C2')
ax3.plot(Saponite, depth, label='Saponite-Fe-Na', color='C3')
ax3.plot(Chamosite, depth, label='Chamosite', color='C9')
#ax3.plot(time_plot, Hydroxyapatite, label='Hydroxyapatite', color='C4')
ax3.plot(Pyrite, depth, label='Pyrite', color='C5')
#ax3.set_xscale('log')
ax3.set_xlabel('Mineral amount [mol/kg]', fontsize=14)
ax3.set_ylabel('Depth [m]', fontsize=14)

labelLines(ax3.get_lines(), zorder=2.5)
#ax3.legend(fontsize=8, loc='upper right')
plt.savefig('Figures/Water_column_t_0.svg', format='svg')
plt.show()
# %%
