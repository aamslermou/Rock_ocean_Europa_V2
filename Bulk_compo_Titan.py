import numpy as np 
import Compo_ocean_Moon as compo
import pandas as pd

# Bulk composition given in kg : 

# 1st method : taking an ocean 700km deep and putting the percentages (not normalized) of cometary composition (cf cahier: to check which one)
# M_CO2 = 1.22e21
# M_CH4 = 5.87e19
# M_NH3 = 3.28e20
# M_Ar  = 2.84e17
# M_Kr  = 2.40e16
# M_Xe  = 1.17e16

# 2nd method : We consider that the ice content of Titan is a % of the total mass. We compute the mass of each volatile as a % of the icy mass (normalized 67P data)
# There're two cases in litt (Lunine 2009) : 30% ice or 50% ice

M_bulk = 1.345e23 # kg                                               # Bulk mass of Titan


M  = pd.DataFrame(np.array([[18.01528e-3, 44.0095e-3, 16.0425e-3, 39.948e-3, 83.798e-3, 131.29e-3, 17.03052e-3, 61.0168e-3, 18.03846e-3, 60.0089e-3, 108.9936e-3]]),      # Molar masses in [kg/mol]

                                    columns= ['H2O', 'CO2', 'CH4', 'Ar', 'Kr', 'Xe', 'NH3', 'HCO3-', 'NH4+', 'CO32-', 'NH2COO-'])

molar_mass = np.array([M['H2O'][0], M['CO2'][0], M['NH3'][0], M['CH4'][0], M['Ar'][0], M['Kr'][0], M['Xe'][0]])
    # 50% ice : 
# M_ice = 0.5 * M_bulk

# M_H2O, M_CO2, M_NH3, M_O2, M_CO, M_CH4, M_N2, M_Ar, M_Kr, M_Xe = compo.compo_67P * M_ice   # 67P composition
# bulk_compo_67P = np.array([M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

# M_CO2, M_NH3, M_CH4, M_Ar, M_Kr, M_Xe = compo.compo_end_CO2rich * M_ice                    # Endmember CO2-rich
# bulk_compo_CO2rich = np.array([M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

# M_CO2, M_NH3, M_CH4, M_Ar, M_Kr, M_Xe = compo.compo_end_CO2poor * M_ice                    # Endmember CO2-poor
# bulk_compo_CO2poor = np.array([M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

    # 30% ice :
M_ice = 0.3 * M_bulk
bulk_mass_67P = np.zeros(7)

print(M_ice)
M_H2O, M_CO2, M_NH3, M_O2, M_CO, M_CH4, M_N2, M_Ar, M_Kr, M_Xe = compo.compo_67P * M_ice   # 67P composition
bulk_compo_67P = np.array([M_H2O,M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])
#print(bulk_compo_67P)

#Total molar mass : 
compo_titan = np.array([compo.compo_67P[0],compo.compo_67P[1],compo.compo_67P[2],compo.compo_67P[5],compo.compo_67P[7],compo.compo_67P[8],compo.compo_67P[9]])
M_tot = sum([molar_mass[i]*compo_titan[i] for i in range(7)])
#print(M_tot)
for i in range(7):
    bulk_mass_67P[i] = bulk_compo_67P[i] * (molar_mass[i]/M_tot)  # 67P composition

V_Ocean = M_H2O / 1000  # m^3
print(V_Ocean)
print(bulk_mass_67P)

M_H2O, M_CO2, M_NH3, M_O2, M_CO, M_CH4, M_N2, M_Ar, M_Kr, M_Xe = compo.compo_end_CO2rich * M_ice             # Endmember CO2-rich
bulk_compo_CO2rich = np.array([M_H2O,M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

M_H2O, M_CO2, M_NH3, M_O2, M_CO, M_CH4, M_N2, M_Ar, M_Kr, M_Xe = compo.compo_end_CO2poor * M_ice             # Endmember CO2-poor
bulk_compo_CO2poor = np.array([M_H2O,M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

bulk_mass_CO2poor = np.zeros(7)
for i in range(7):
    bulk_mass_CO2poor[i] = bulk_compo_CO2poor[i] * (molar_mass[i]/molar_mass[0]) 
#print(bulk_mass_CO2poor)
    # 40% ice :
# M_ice = 0.4 * M_bulk

# M_H2O, M_CO2, M_NH3, M_O2, M_CO, M_CH4, M_N2, M_Ar, M_Kr, M_Xe = compo.compo_67P * M_ice           # 67P composition
# bulk_compo_67P = np.array([M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

# M_H2O, M_CO2, M_NH3, M_O2, M_CO, M_CH4, M_N2, M_Ar, M_Kr, M_Xe = compo.compo_end_CO2rich * M_ice   # Endmember CO2-rich
# bulk_compo_CO2rich = np.array([M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

# M_H2O, M_CO2, M_NH3, M_O2, M_CO, M_CH4, M_N2, M_Ar, M_Kr, M_Xe = compo.compo_end_CO2poor * M_ice   # Endmember CO2-poor
# bulk_compo_CO2poor = np.array([M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])

#print(bulk_compo_CO2poor)

# Bulk compo with Enceladus data
M_H2O, M_CO2, M_NH3, M_CH4, M_Ar, M_Kr, M_Xe = compo.compo_Enceladus_maxNH3 * M_ice   # Enceladus grains composition
bulk_compo_Enceladus_maxNH3 = np.array([M_H2O,M_CO2,M_NH3,M_CH4,M_Ar,M_Kr,M_Xe])
#print(bulk_compo_67P)

#Total molar mass : 
M_tot = sum([molar_mass[i]*compo.compo_Enceladus_maxNH3[i] for i in range(7)])
#print(M_tot)

bulk_mass_Enceladus_maxNH3 = np.zeros(7)

for i in range(7):
    bulk_mass_Enceladus_maxNH3[i] = bulk_compo_Enceladus_maxNH3[i] * (molar_mass[i]/M_tot)  

V_Ocean = M_H2O / 1000  # m^3

