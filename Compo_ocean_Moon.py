# This file gathers the initial compositions of Europa's primordial ocean for different formation scenarios : 
# Order of the molecules : H2O, CO2, NH3, O2, CO, CH4, N2, Ar, Kr, Xe
import numpy as np 
import pandas as pd
# COMPO BASED ON 67P from Rubin et al. 2019 :

compo_67P = np.array([1, 0.047, 0.0067, 0.031, 0.031, 0.0034, 0.00089, 0.0000058, 0.00000049, 0.00000024]) #Not normalized 
compo_67P = compo_67P/compo_67P.sum()
compo_67P = pd.DataFrame(compo_67P.reshape(1,10),
                         
                         columns= ['H2O', 'CO2', 'NH3', 'O2', 'CO', 'CH4', 'N2', 'Ar', 'Kr', 'Xe'])
#print(compo_67P)

# Endmember CO2-rich composition from Bockelée-Morvan 2017 : 

compo_end_CO2rich = np.array([1, 0.3, 0.0067, 0.031, 0.031, 0.0034, 0.00089, 0.0000058, 0.00000049, 0.00000024])
compo_end_CO2rich = compo_end_CO2rich/compo_end_CO2rich.sum()
#print(compo_end)

# Endmember CO2-rich composition from Bockelée-Morvan 2017 : 

compo_end_CO2poor = np.array([1, 0.025, 0.0067, 0.031, 0.031, 0.0034, 0.00089, 0.0000058, 0.00000049, 0.00000024])
compo_end_CO2poor = compo_end_CO2poor/compo_end_CO2poor.sum()

# Composition based on excel file from Kathy :

compo_Kathy_1 = [0.941855074, 0.025, 0.0015, 0.025479746, 0.006, 0.00016, 4.60*1e-6, 3.90*1e-7, 1.90*1e-7]
compo_Kathy_2 = [0.586294015, 0.369, 0.0006, 0.0363380805, 0.007, 0.00072, 4.60*1e-6, 3.90*1e-7, 1.90*1e-7]
compo_Kathy_3 = [0.839187155, 0.14, 0.0006, 0.0171187665, 0.0023, 0.00072, 4.60*1e-6, 3.90*1e-7, 1.90*1e-7]
compo_Kathy_4 = [0.641152666, 0.32, 0.0006, 0.032822154, 0.0047, 0.00072, 4.60*1e-6, 3.90*1e-7, 1.90*1e-7]

compo_Kathy = np.array([compo_Kathy_1,
                        compo_Kathy_2,
                        compo_Kathy_3,
                        compo_Kathy_4]) 

compo_test = [0.947356074, 0.025, 0.033, 0.0000008, 0.00006, 1.6e-9, 4.60*1e-9, 3.90*1e-7, 1.90*1e-7]

# Compo Enceladus measured by NIMS on board Cassini 
# H2O, CO2, NH3, CH4, Ar, Kr, Xe
compo_Enceladus_minNH3 = np.array([0.985, 0.008, 0.004, 0.003,  0.0000058, 0.00000049, 0.00000024])
compo_Enceladus_minH3 = compo_Enceladus_minNH3/compo_Enceladus_minNH3.sum()

compo_Enceladus_maxNH3 = np.array([0.976, 0.008, 0.013, 0.003,  0.0000058, 0.00000049, 0.00000024])
compo_Enceladus_maxNH3 = compo_Enceladus_maxNH3/compo_Enceladus_maxNH3.sum()
