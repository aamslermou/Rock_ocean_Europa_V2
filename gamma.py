#%%                                     ########################### Import packages ##########################
import numpy as np
import pandas as pd

from search_function import getIndexes_onecolumn
#%%                                         ############################ DATA #################################

                    
species = pd.DataFrame(np.array([['H2O', 647.3, 22.12, 0.3434, 18.01528, 0, 0.92, 1.40],
                                 ['CO2', 304.19, 7.383, 0.224,  44.0095, 0, 0.75, 2.45],
                                 ['NH3', 405.6, 11.28, 0.25, 17.03052, 0, 1.6292, 2.9852],
                                 ['HCO3-', None , None , None , 61.0168, -1, 8.0756, 8.6806],
                                 ['H+', None, None, None, 1.008, 1, 0.1378, 0.1e-15],
                                 ['NH4+', None, None, None, 18.03846, 1, 4.8154, 4.6028],
                                 ['CO32-',None, None, None, 60.0089, -2, 10.828, 10.769],
                                 ['OH-',None, None, None, 17.00734, -1, 1.1, 1],
                                 ['NH2COO-', None, None, None, 108.9936, -1, 4.3022, 4.1348]]),

                       
                        columns = ['formula', 'Tc', 'Pc', 'w', 'M', 'charge', 'r_coeff', 
                           'q_coeff'])

# u0_coeff = pd.DataFrame(np.array([[0, 8.8383,594.72,577.05,52.7305,361.39,28.2779],
#                                  [8.8383,302.25,2500.0,526.305,-424.01,2500.0,2500.0],
#                                  [594.72,2500,1090.8,534.01,785.98,524.13,498.15],
#                                  [577.05,526.305,534.01,771.04,505.55,800.01,613.25],
#                                  [52.7305,-424.01,785.98,505.55,0,226.60,44.849],
#                                  [361.39,2500.0,524.13,800.01,226.6,1458.3,2500.0],
#                                  [28.2779, 2500.0, 498.15, 613.25, 44.849, 2500.0, 3343.1]]),
                       
#                         columns = ['H2O', 'CO2', 'NH3', 'HCO3-','NH4+','CO32-','NH2COO-'])
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

# ut_coeff = pd.DataFrame(np.array([[0,0.86293,7.1827,-0.38795,0.50922,3.3516,8.0238],
#                                  [0.8629,0.35870,0,-3.7340,8.6951,0,0],
#                                  [7.1827,0,7.0912,5.3111,6.1271,4.9305,6.6532],
#                                  [-0.388,-3.734,5.3111,-0.01981,-0.00795,1.7241,3.0580],
#                                  [0.5092,8.6951,6.1227,-0.00795,0,4.0555,12.047],
#                                  [3.3516,0,4.9305,1.724,4.056,-1.34,0],
#                                  [8.0238,0,6.6532,3.0580,12.047,0,-15.920]]),
                       
#                         columns = ['H2O', 'CO2', 'NH3','HCO3-','NH4+','CO32-','NH2COO-'])
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
N = 9
psi = np.zeros((N,N))
u   = np.zeros((N,N))
T   = 333
    
for i in range(0,N):
   for j in range(0,N):
      u[i,j]   = u0_coeff.iat[i,j] + ut_coeff.iat[i,j]*(T-298.15)   # u computed with data from litterature (darde et al. 2010)

for i in range(0,N):
   for j in range(0,N):
      psi[i,j] = np.exp(-(u[i,j]-u[j,j])/(T))    

#print(psi)
#%%                        #################### Functions to compute the activity coefficient of water ########################

                        
from search_function import getIndexes_onecolumn

def gamma_solvent_uniquac(T, x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm, component, species, psi):

    
# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# psi        : psi constants computed from u0 and ut coefficient according to UNIQUAC model
 
# OUTPUTS : 
# gamma_UNIQUAC_water : Coefficient of activity of water computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity for water, using the extended UNIQUAC model presented
# in Marounina's thesis p34-36. 

                                #####################################################
    
# Extraction of constants : 
    x       = [x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm]
    N       = len(x)
    z       = 10                    # Coordination number. Fixed as 10 according to Maurer and Prausnitz (1978)
    index   = getIndexes_onecolumn(species,component)[0]
    r_coeff = species['r_coeff']    # Molecular volume parameter (Darde et al. 2010)
    q_coeff = species['q_coeff']    # Molecular surface parameter (Darde et al. 2010)

# Computation of necessary constants :

    phi       = []
    theta     = []
    sum_phi   = 0
    sum_theta = 0

    for i in range(0,N):
        sum_phi   += x[i]*r_coeff[i]
        sum_theta += x[i]*q_coeff[i]

    for i in range(0,N):
        phi.append(x[i]*r_coeff[i]/sum_phi)
        theta.append(x[i]*q_coeff[i]/sum_theta)
    
    # Constants specific to water :     

    phi_i   = phi[index]
    theta_i = theta[index]
    x_i     = x[index]
    q_i     = q_coeff[index]

# Computation of symmetric coefficients of activity : 
 
    log_gamma_c_w = np.log(phi_i/x_i) + 1 - phi_i/x_i - 5*q_i*(np.log(phi_i/theta_i) + 1 - phi_i/theta_i)

    sum_mix       = []
    sum_tot       = 0 
    sum_theta_psi = 0
    pre_sum_mix   = 0

    for j in range(0,N):     
        for k in range(0,N):
            pre_sum_mix += theta[k]*psi[k,j]
        sum_mix.append(pre_sum_mix)

    for j in range(0,N):
        sum_tot += theta[j]*psi[index,j]/sum_mix[j]
        sum_theta_psi += theta[j]*psi[j,index]
    
    log_gamma_r_w = q_i*(1 - np.log(sum_theta_psi) - sum_tot)

# Computation of Debye_Huckel term 

    A      = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2    # VALIDE BETWEEN [273.15-373.15] K
    M      = species['M']
    M_i    = M[index]            # Molar mass of considered specie
    b      = 1.5                 # Constant [kg**1/2 mol**-1/2]
    z_dh   = species['charge']

    sum_ions    = 0
    sum_solvent = 0

    for i in range(0,N):
        sum_ions += x[i]*(z_dh[i]**2)
        if z_dh[i]==0:
            sum_solvent += x[i]*(M[i]*1e-3)
    

    # I = 0.5*sum_ions/sum_solvent   # Ionic strength
    # print(I) 
    I = 0.1
    X = b*I**(1/2)
    log_gamma_DH_w = (2/3)*(M_i*1e-3)*A*I**(3/2)*((3/X**3)*(1 + X - 1/(1 + X) -2*np.log(1 + X)))

# Final gamma computation

    # print(log_gamma_c_w,log_gamma_r_w,log_gamma_DH_w)

    gamma_UNIQUAC_w = np.exp(log_gamma_c_w + log_gamma_r_w + log_gamma_DH_w) 

    return gamma_UNIQUAC_w


#%%                                     ################################ Test ###############################
# x0 = [0.9, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
# T = 300
# P = 1e5

# res = gamma_solvent_uniquac(T,x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7], x0[8],'H2O',species,psi)
# print(res)
# %%                       #################### Functions to compute the activity coefficient of solutes ########################

 
def gamma_UNIQUAC(T, x_H2O, x_CO2, x_NH3, x_HCO3m, x_NH4p, x_CO32m, x_NH2COOm, species, psi):
    
# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# psi        : psi constants computed from u0 and ut coefficient according to UNIQUAC model

# OUTPUTS : 
# gamma_UNIQUAC : Coefficient of activity of any solute or ionic species computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity for water, using the extended UNIQUAC model presented
# in Marounina's thesis p34-36. 

                                #####################################################
    
# Extraction of constants :   
    x           = [x_H2O,x_CO2,x_NH3,x_HCO3m,x_NH4p,x_CO32m,x_NH2COOm]

    cmps = range(len(x))

    N           = len(x)
    z           = 10        # Coordination number. Fixed as 10 according to Maurer and Prausnitz (1978)
    r_coeff     = species['r_coeff']
    q_coeff     = species['q_coeff']
    index_water = getIndexes_onecolumn(species,'H2O')[0]

    taus = psi

# Computation of necessary constants :

    rsxs   = [r_coeff[i]*x[i] for i in cmps]
    rsxs_sum_inv = 1.0/sum(rsxs)
    phis   = [rsxs[i]*rsxs_sum_inv for i in cmps]

    qsxs   = [q_coeff[i]*x[i] for i in cmps]
    qsxs_sum_inv = 1.0/sum(qsxs)
    thetas = [qsxs[i]*qsxs_sum_inv for i in cmps]

    Ss     = [sum([thetas[j]*taus[j][i] for j in cmps]) for i in cmps]
    VsSs   = [thetas[j]/Ss[j] for j in cmps]

#Computation of symmetric coefficients of activity : 

    loggammac = []
    loggammar = []
    for i in cmps:
        x1 = phis[i]/x[i]
        x2 = phis[i]/thetas[i]
        loggammac.append(np.log(x1) + 1.0 - x1 - 5.0*q_coeff[i]*(np.log(x2) + 1.0 - x2))
        loggammar.append(q_coeff[i]*(1.0 - np.log(Ss[i]) - sum([taus[i][j]*VsSs[j] for j in cmps])))
    
    q_w  = q_coeff[index_water]
    r_w  = r_coeff[index_water]  
    
# Computation of infinite dilution terms 

    loggammacinf = []
    loggammarinf = []
    for i in cmps:
        loggammacinf.append(np.log(r_coeff[i]/r_w) + 1 - (r_coeff[i]/r_w) - 0.5*z*q_coeff[i]*(np.log(r_coeff[i]*q_w/(r_w*q_coeff[i])) + 1 - r_coeff[i]*q_w/(r_w*q_coeff[i])))
        loggammarinf.append(q_coeff[i]*(1 - np.log(taus[index_water,i]) - taus[i,index_water]))

# Computation of Debye_Huckel term 

    A     = 1.31 + 1.335*1e-3*(T - 273.15) + 1.164*1e-5*(T - 273.15)**2    # VALIDE BETWEEN [273.15-373.15] K
    Mw    = 18*1e-3                 # Molar masse of water
    b     = 1.5                     # Constant [kg**1/2 mol**-1/2]
    z_dh  = species['charge']
    
    m = np.zeros(N)
    I = 0

    for i in range(N):
        m[i] = x[i]/(Mw*x_H2O)
        I += 0.5*m[i]*z_dh[i]**2    # Ionic strength

    loggammaDH = []
    for i in cmps:
        loggammaDH.append(-(z_dh[i]**2)*(A*I**0.5)/(1 + b*I**0.5))

# Final gamma computation

    gamma_UNIQUAC = []
    for i in cmps:
        gamma_UNIQUAC.append(np.exp(loggammac[i] - loggammacinf[i] + loggammar[i] - loggammarinf[i] + loggammaDH[i])) 

    return gamma_UNIQUAC                                                       
#%%                                     ################################ Test ###############################

# x0 = [0.96,1e-4, 6*1e-3 , 9*1e-3, 2.271*1e-2, 3.6*1e-3, 6.4*1e-3]
# T = 333
# P = 1e5

# res = gamma_UNIQUAC(T,x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6],species,psi)

# print(res)
 # %%                                       ##################### Gamma 2 species water ####################
from search_function import getIndexes_onecolumn

def gamma_water_uniquac_2species(T, x_H2O, x_CO2, component, species, psi):

# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# psi        : psi constants computed from u0 and ut coefficient according to UNIQUAC model
 
# OUTPUTS : 
# gamma_UNIQUAC_water : Coefficient of activity of water computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity for water, using the extended UNIQUAC model presented
# in Marounina's thesis p34-36. 

                                #####################################################
    
# Extraction of constants : 
    x       = [x_H2O, x_CO2]
    N       = len(x)
    z       = 10                    # Coordination number. Fixed as 10 according to Maurer and Prausnitz (1978)
    index   = getIndexes_onecolumn(species,component)[0]
    r_coeff = species['r_coeff']    #Molecular volume parameter (Darde et al. 2010)
    q_coeff = species['q_coeff']    #Molecular surface parameter (Darde et al. 2010)


# Computation of necessary constants :

    phi       = []
    theta     = []
    sum_phi   = 0
    sum_theta = 0

    for i in range(0,N):
        sum_phi   += x[i]*r_coeff[i]
        sum_theta += x[i]*q_coeff[i]

    for i in range(0,N):
        phi.append(x[i]*r_coeff[i]/sum_phi)
        theta.append(x[i]*q_coeff[i]/sum_theta)
    
    # Constants specific to water :     

    phi_i   = phi[index]
    theta_i = theta[index]
    x_i     = x[index]
    q_i     = q_coeff[index]

# Computation of symmetric coefficients of activity : 
 
    log_gamma_c_w = np.log(phi_i/x_i) + 1 - phi_i/x_i - 5*q_i*(np.log(phi_i/theta_i) + 1 - phi_i/theta_i)

    sum_mix       = []
    sum_tot       = 0 
    sum_theta_psi = 0
    pre_sum_mix   = 0

    for j in range(0,N):     
        for k in range(0,N):
            pre_sum_mix += theta[k]*psi[k,j]
        sum_mix.append(pre_sum_mix)

    for j in range(0,N):
        sum_tot += theta[j]*psi[index,j]/sum_mix[j]
        sum_theta_psi += theta[j]*psi[j,index]
    
    log_gamma_r_w = q_i*(1 - np.log(sum_theta_psi) - sum_tot)


# Computation of Debye_Huckel term 

    # A   = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2    # VALIDE BETWEEN [273.15-373.15] K
    # M_w = 18.01528            # Molar mass of water
    # b   = 1.5                 # Constant [kg**1/2 mol**-1/2]
    # z   = species['charge']
    # m   = []
    # for i in range(0,N):
    #     m.append(x[i]/(M_w*x_H2O))      # molalités : rajouter une input avec les concentrations ? Relier à xi ?? 
    # I   = 0.5*sum(m*(z**2))   # Ionic strength
    # #print(I)
    # X = b*I**(1/2)
    # log_gamma_DH_w = (2/3)*M_w*A*I**(3/2)*((3/X**3)*(1 + X - 1/(1 + X) -2*np.log(1 + X)))


# Final gamma computation

    #print(np.exp(log_gamma_c_w),np.exp(log_gamma_r_w),np.exp(log_gamma_DH_w))

    gamma_UNIQUAC_w = np.exp(log_gamma_c_w + log_gamma_r_w) 

    return gamma_UNIQUAC_w

# %%                                      ################## Gamma computation for 2 species #################
def gamma_UNIQUAC_2species(T, x_H2O, x_CO2, component, species, psi):
    
# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!


# psi        : psi constants computed from u0 and ut coefficient according to UNIQUAC model

# OUTPUTS : 
# gamma_UNIQUAC : Coefficient of activity of any solute or ionic species computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity for water, using the extended UNIQUAC model presented
# in Marounina's thesis p34-36. 

                                #####################################################
    
# Extraction of constants :   
    x           = [x_H2O,x_CO2]
    N           = len(x)
    z           = 10        # Coordination number. Fixed as 10 according to Maurer and Prausnitz (1978)
    r_coeff     = species['r_coeff']
    q_coeff     = species['q_coeff']
    index       = getIndexes_onecolumn(species,component)[0]
    index_water = getIndexes_onecolumn(species,'H2O')[0]

# Computation of necessary constants :

    phi       = []
    theta     = []
    sum_phi   = 0
    sum_theta = 0

    for i in range(0,N):
        sum_phi   += x[i]*r_coeff[i]
        sum_theta += x[i]*q_coeff[i]

    for i in range(0,N):
        phi.append(x[i]*r_coeff[i]/sum_phi)
        theta.append(x[i]*q_coeff[i]/sum_theta)

    phi_i   = phi[index]
    theta_i = theta[index]
    x_i     = x[index]
    q_i     = q_coeff[index]
    r_i     = r_coeff[index]
    q_w     = q_coeff[index_water]
    r_w     = r_coeff[index_water]  
    psi_wi  = psi[index_water,index]
    psi_iw  = psi[index,index_water]
    
    #print(phi_i,theta_i,x_i,q_i,q_w,psi_wi,psi_iw)

        
# Computation of symmetric coefficients of activity : 

    log_gamma_c = np.log(phi_i/x_i) + 1 - phi_i/x_i - 0.5*z*q_i*(np.log(phi_i/theta_i) + 1 - phi_i/theta_i)

    sum_mix = []
    sum_tot = 0
    sum_theta_psi = 0
    pre_sum_mix = 0 


    for j in range(0,N):     
        for k in range(0,N):
            pre_sum_mix += theta[k]*psi[k,j]
        sum_mix.append(pre_sum_mix)

    for j in range(0,N):
        sum_tot += theta[j]*psi[index,j]/sum_mix[j]
        sum_theta_psi += theta[j]*psi[j,index]
    
    log_gamma_r   = q_i*(1 - np.log(sum_theta_psi) - sum_tot)

# Computation of infinite dilution terms 

    log_gamma_c_inf = np.log(r_i/r_w) + 1 - (r_i/r_w) - 0.5*z*q_i*(np.log(r_i*q_w/(r_w*q_i)) + 1 - r_i*q_w/(r_w*q_i))

    log_gamma_r_inf = q_i*(1 - np.log(psi_wi) - psi_iw)

    # print(np.exp(log_gamma_c_inf),np.exp(log_gamma_r_inf))
# Computation of Debye_Huckel term 

        # A   = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2    # VALIDE BETWEEN [273.15-373.15] K
        # M_w = 18.01528            # Molar mass of water
        # b   = 1.5                 # Constant [kg**1/2 mol**-1/2]
        # z   = species['charge']
        # m   = []
        # for i in range(0,N):
        #     m.append(x[i]/(M_w*x_H2O))
        # I   = 0.5*sum(m*(z**2))   # Ionic strength
        # #print(I)
        # log_gamma_DH = -z[index]**2*(A*I**0.5)/(1 + b*I**0.5)


# Final gamma computation

    gamma_UNIQUAC = np.exp(log_gamma_c - log_gamma_c_inf + log_gamma_r - log_gamma_r_inf) 

    return gamma_UNIQUAC  


# %%                                              ################### New gamma #################

from thermo import UNIQUAC

def gamma_new(T,x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm,species,u):

# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# u          : u constants computed from u0 and ut coefficient according to UNIQUAC model
 
# OUTPUTS : 
# gamma : Coefficient of activity computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity using the thermo package, adding the Debye-Hückel term

    x = [x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm]
    N = len(x)
    
    Aij0 = np.zeros((N,N))

    for i in range(0,N):
        for j in range(0,N):       
            Aij0[i,j] = -(u[i,j]-u[j,j])    

# Definition of constants for UNIQUAC

    tausA = tausC = tausD = tausE = tausF = [[0.0]*N for i in range(N)]

    tausB = Aij0

    ABCDEF = (tausA, tausB, tausC, tausD, tausE,  tausF)

    rs = species['r_coeff']
    qs = species['q_coeff']

    GE = UNIQUAC(T=T, xs=x, rs=rs, qs=qs, ABCDEF=ABCDEF)

    gamma = GE.gammas()

    Z   = species['charge']
    I   = 0 
    M_w = 18.01528*1e-3     
    A   = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2    # VALIDE BETWEEN [273.15-373.15] K
    b   = 1.5 

    # for i in range(N):
    #     mi = x[i]/(M_w*x_H2O)
    #     I += 0.5*mi*Z[i]**2
        
    # for i in range(N):
    #     if i == 0:
    #         X = b*I**(0.5)
    #         gamma[0] = np.exp(np.log(gamma[0]) + (2/3)*M_w*A*I**(3/2)*((3/X**3)*(1 + X - 1/(1 + X) -2*np.log(1 + X))))
    #     else:
    #         gamma[i] = np.exp(np.log(gamma[i]) + -(Z[i]**2)*(A*I**(0.5))/(1 + b*I**(0.5)))

    return gamma


# %%                                              ################# Test new gamma ###############

# gamma = gamma_new(T,0.9, 0.05, 1e-10, 1e-10, 1e-10, 0.05, 1e-10, 1e-10, 1e-10,species,u)
# print(gamma)

# %%                                                  ################# Plots ###################

# T = 333

# gamma_H2O_list     = []
# gamma_CO2_list     = []
# gamma_NH3_list     = []
# gamma_NH4p_list    = []
# gamma_Hp_list      = []
# gamma_HCO3m_list   = []
# gamma_CO32m_list   = []
# gamma_NH2COOm_list = []
# gamma_OHm_list     = []
# gamma_list_H2O     = np.zeros((100,9))
# gamma_list_CO2     = np.zeros((100,9))
# gamma_list_NH3     = np.zeros((100,9))
# gamma_list_HCO3m   = np.zeros((100,9))
# gamma_list_Hp      = np.zeros((100,9))
# gamma_list_NH4p    = np.zeros((100,9))
# gamma_list_OHm     = np.zeros((100,9))
# gamma_list_CO32m   = np.zeros((100,9))
# gamma_list_NH2COOm = np.zeros((100,9))

# x_list = np.linspace(1e-5,0.98,100)
# i = 0
# #for i in range(0,10):
# for x in x_list:

#     gamma_list_H2O[i,:]     = (gamma_new(T, x, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, species, u))
#     gamma_list_CO2[i,:]     = (gamma_new(T, (1-x)/8, x, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, species, u))
#     gamma_list_NH3[i,:]     = (gamma_new(T, (1-x)/8, (1-x)/8, x, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, species, u))
#     gamma_list_HCO3m[i,:]   = (gamma_new(T, (1-x)/8, (1-x)/8, (1-x)/8, x, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, species, u))
#     gamma_list_Hp[i,:]      = (gamma_new(T, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, x, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, species, u))
#     gamma_list_NH4p[i,:]    = (gamma_new(T, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, x, (1-x)/8, (1-x)/8, (1-x)/8, species, u))
#     gamma_list_OHm[i,:]     = (gamma_new(T, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, x, (1-x)/8, (1-x)/8, species, u))
#     gamma_list_CO32m[i,:]   = (gamma_new(T, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, x, (1-x)/8, species, u))
#     gamma_list_NH2COOm[i,:] = (gamma_new(T, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, (1-x)/8, x, species, u))

#     i += 1

# plt.figure(1)
# gamma_H2O_list     = gamma_list_H2O[:,0]
# gamma_CO2_list     = gamma_list_CO2[:,1]
# gamma_NH3_list     = gamma_list_NH3[:,2]
# gamma_HCO3m_list   = gamma_list_HCO3m[:,3]
# gamma_Hp_list      = gamma_list_Hp[:,4]
# gamma_NH4p_list    = gamma_list_NH4p[:,5]
# gamma_OHm_list     = gamma_list_OHm[:,6]
# gamma_CO32m_list   = gamma_list_CO32m[:,7]
# gamma_NH2COOm_list = gamma_list_NH2COOm[:,8]

# #plt.plot(x_list,gamma_H2O_list, color='red', label='$\gamma_{H2O}$')
# plt.plot(x_list,gamma_CO2_list, color='green', label='$\gamma_{CO2}$')
# plt.plot(x_list,gamma_NH3_list, color='blue', label='$\gamma_{NH3}$')
# plt.plot(x_list,gamma_HCO3m_list, color='purple', label='$\gamma_{HCO3^-}$')
# plt.plot(x_list,gamma_Hp_list, color='pink', label='$\gamma_{H^+}$')
# plt.plot(x_list,gamma_NH4p_list, color='yellow', label='$\gamma_{NH4^+}$')
# plt.plot(x_list,gamma_CO32m_list, color='black', label='$\gamma_{CO32^-}$')
# plt.plot(x_list,gamma_NH2COOm_list, color='orange', label='$\gamma_{NH2COO^-}$')
# plt.plot(x_list,gamma_OHm_list, color='cyan', label='$\gamma_{OH^-}$')

# plt.title('T = 333 K')
# plt.ylabel('$\gamma$')
# plt.legend()

# plt.show()
# %%                                               ########## Gamma for reduced systems ###########""
from thermo import UNIQUAC

def gamma_new_3species(T,x_H2O, x_CO2, x_NH3,species,u):

# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# u          : u constants computed from u0 and ut coefficient according to UNIQUAC model
 
# OUTPUTS : 
# gamma : Coefficient of activity computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity using the thermo package, adding the Debye-Hückel term

    x = [x_H2O, x_CO2, x_NH3]
    N = len(x)
    
    Aij0 = np.zeros((N,N))

    for i in range(0,N):
        for j in range(0,N):       
            Aij0[i,j] = -(u[i,j]-u[j,j])    

# Definition of constants for UNIQUAC

    tausA = tausC = tausD = tausE = tausF = [[0.0]*N for i in range(N)]

    tausB = Aij0

    ABCDEF = (tausA, tausB, tausC, tausD, tausE,  tausF)

    species3 = species[0:3]
    rs = species3['r_coeff']
    qs = species3['q_coeff']

    GE = UNIQUAC(T=T, xs=x, rs=rs, qs=qs, ABCDEF=ABCDEF)

    gamma = GE.gammas()

    return gamma


def gamma_new_8species(T,x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_OHm, x_NH2COOm, species,u):

# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# u          : u constants computed from u0 and ut coefficient according to UNIQUAC model
 
# OUTPUTS : 
# gamma : Coefficient of activity computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity using the thermo package, adding the Debye-Hückel term

    x = [x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_OHm, x_NH2COOm]
    N = len(x)

    Aij0 = np.zeros((N,N))

    for i in range(0,N):
        for j in range(0,N):       
            Aij0[i,j] = -(u[i,j]-u[j,j])    

# Definition of constants for UNIQUAC

    tausA = tausC = tausD = tausE = tausF = [[0.0]*N for i in range(N)]

    tausB = Aij0

    ABCDEF = (tausA, tausB, tausC, tausD, tausE,  tausF)

    rs = species['r_coeff']
    qs = species['q_coeff']

    GE = UNIQUAC(T=T, xs=x, rs=rs, qs=qs, ABCDEF=ABCDEF)

    gamma = GE.gammas()

    Z   = species['charge']
    I   = 0 
    M_w = 18.01528*1e-3     
    A   = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2    # VALIDE BETWEEN [273.15-373.15] K
    b   = 1.5 

    for i in range(N):
        mi = x[i]/(M_w*x_H2O)
        I += 0.5*mi*Z[i]**2
    
    for i in range(N):
        if i == 0:
            X = b*I**(0.5)
            gamma[0] = np.exp(np.log(gamma[0]) + (2/3)*M_w*A*I**(3/2)*((3/X**3)*(1 + X - 1/(1 + X) -2*np.log(1 + X))))
        else:
            gamma[i] = np.exp(np.log(gamma[i]) + -(Z[i]**2)*(A*I**(0.5))/(1 + b*I**(0.5)))
    return gamma


def gamma_new_5species(T,x_H2O, x_CO2, x_NH3, x_Hp, x_OHm, species,u):

# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# u          : u constants computed from u0 and ut coefficient according to UNIQUAC model
 
# OUTPUTS : 
# gamma : Coefficient of activity computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity using the thermo package, adding the Debye-Hückel term

    x = [x_H2O, x_CO2, x_NH3, x_Hp, x_OHm]
    N = len(x)
    
    Aij0 = np.zeros((N,N))

    for i in range(0,N):
        for j in range(0,N):       
            Aij0[i,j] = -(u[i,j]-u[j,j])    

# Definition of constants for UNIQUAC

    tausA = tausC = tausD = tausE = tausF = [[0.0]*N for i in range(N)]

    tausB = Aij0

    ABCDEF = (tausA, tausB, tausC, tausD, tausE,  tausF)

    rs = species['r_coeff']
    qs = species['q_coeff']

    GE = UNIQUAC(T=T, xs=x, rs=rs, qs=qs, ABCDEF=ABCDEF)

    gamma = GE.gammas()

    Z   = species['charge']
    I   = 0 
    M_w = 18.01528*1e-3     
    A   = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2    # VALIDE BETWEEN [273.15-373.15] K
    b   = 1.5 

    for i in range(N):
        mi = x[i]/(M_w*x_H2O)
        I += 0.5*mi*Z[i]**2
    
    for i in range(N):
        if i == 0:
            X = b*I**(0.5)
            gamma[0] = np.exp(np.log(gamma[0]) + (2/3)*M_w*A*I**(3/2)*((3/X**3)*(1 + X - 1/(1 + X) -2*np.log(1 + X))))
        else:
            gamma[i] = np.exp(np.log(gamma[i]) + -(Z[i]**2)*(A*I**(0.5))/(1 + b*I**(0.5)))
    return gamma

# %%                                                   ############## Gamma Krop ###############
def gamma_Krop(T,x_H2O,x_CO2,x_NH3,x_HCO3m,x_NH4p,x_CO32m,x_NH2COOm,species):

    Z = species['charge']
    A   = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2
    I   = 0 
    M_w = 18.01528*1e-3    
  
    x = [x_H2O, x_CO2, x_NH3, x_HCO3m, x_NH4p, x_CO32m, x_NH2COOm]  
    N = len(x)

    m = np.zeros(N)
    for i in range(N):
        m[i] = x[i]/(M_w*x_H2O)
        I += 0.5*m[i]*Z[i]**2      

    
    m_H2O, m_CO2, m_NH3, m_HCO3m, m_NH4p, m_CO32m, m_NH2COOm = m

    f1_g = -A*(I**0.5/(1 + 1.2*I**0.5) + 2/1.2*np.log(1 + 1.2*I**0.5))
    f2_g = 1/(4*I**2)*(1 - (1 + 2*I**0.5 + 2*I)*np.exp(-2*I**0.5))

    f2   = 1/(2*I)*(1 - (1 + 2*I**0.5)*np.exp(-2*I**0.5))

    f1_phi = -A*2*I**1.5/(1+1.2*I**0.5)
    f2_phi = np.exp(-2*I**0.5)

    B0_AA = 0.04166 - 9.6975*1e-5*T
    B0_AC = B0_AA 
    B0_CA = -B0_AA*2.2426
    B1_AC = 0.66838 + 6.081*1e-4*T
    B1_CA = -B1_AC*0.18488
    B2_AC = 7.9545*1e-4
    B2_CA = -B2_AC

    m_A = m_NH3 + m_NH4p + m_NH2COOm
    m_C = m_CO2 + m_HCO3m + m_CO32m + m_NH2COOm

    gamma_CO2     = np.exp(2*m_A*(B0_AC + f2*B1_AC) + 3*(m_A**2*B2_AC + 2*m_A*m_C*B2_CA))
    gamma_NH3     = np.exp(2*(m_A*B0_AA + m_C*(B0_CA + B1_CA*f2)) + 3*(2*m_A*m_C*B2_AC + m_C**2*B2_CA))
    gamma_HCO3m   = np.exp(f1_g + 2*(m_A*B0_AA + m_A*(B0_CA + f2*B1_AC) - f2_g*(m_A*m_C*B1_AC + m_A*m_C*B1_CA) + 3*(m_A**2*B2_AC + 2*m_A*m_C*B2_CA)))
    gamma_Hp      = np.exp(f1_g - f2_g*(m_A*m_C*B1_AC + m_C*m_A*B1_CA))
    gamma_OHm     = gamma_Hp
    gamma_NH4p    = np.exp(f1_g + 2*(m_A*B0_AA + m_A*(B0_CA + B1_CA*f2) + m_C*(B0_CA + B1_CA*f2)) - f2_g*(m_A*m_C*B1_AC + m_C*m_A*B1_CA) + 3*(m_C**2*B2_CA + 2*m_A*m_C*B2_AC))
    gamma_CO32m   = np.exp(4*f1_g + 2*m_A*(B0_AC + B1_AC*f2) - 4*f2_g*(m_A*m_C*B1_AC + m_C*m_A*B1_CA) + 3*(m_A**2*B2_AC + 2*m_A*m_C*B2_CA))
    gamma_NH2COOm = np.exp(f1_g + 2*(m_A*B0_AA + m_A*(B0_AC + f2*B1_AC) + m_C*(B0_CA + f2*B1_CA) - f2_g*(m_A*m_C*B1_AC + m_C*m_A*B1_CA)) + 3*(m_A**2*B2_AC + 2*m_A*m_C*(B2_AC + B2_CA) + m_C**2*B2_CA))
    
    a_H2O         = np.exp(-M_w*(f1_phi + m_A**2*B0_AA + m_A*m_C*(B0_CA + B0_AC + (B1_AC + B1_CA)*f2_phi) + 6*m_A*m_C*(m_A*B2_AC + m_C*B2_CA)) - np.log(1 + M_w*sum(m)))
    gamma_H2O     = a_H2O/x_H2O

    gamma = [gamma_H2O, gamma_CO2 ,gamma_NH3, gamma_HCO3m, gamma_Hp, gamma_NH4p, gamma_CO32m, gamma_OHm, gamma_NH2COOm]

    return gamma

def gamma_Krop_3species(T,x_H2O,x_CO2,x_NH3,species):

    Z = species['charge']
    A   = 1.31 + 1.335*10**(-3)*(T - 273.15) + 1.164*10**(-5)*(T - 273.15)**2
    I   = 0 
    M_w = 18.01528*1e-3    
  
    x = [x_H2O, x_CO2, x_NH3]  
    N = len(x)

    m = np.zeros(N)
    for i in range(N):
        m[i] = x[i]/(M_w*x_H2O) 

    m_H2O, m_CO2, m_NH3 = m

    # f1_g = -A*(I**0.5/(1 + 1.2*I**0.5) + 2/1.2*np.log(1 + 1.2*I**0.5))
    # f2_g = 1/(4*I**2)*(1 - (1 + 2*I**0.5 + 2*I)*np.exp(-2*I**0.5))

    # f2   = 1/(2*I)*(1 - (1 + 2*I)*np.exp(-2*I**0.5))

    # f1_phi = -A*2*I**1.5/(1+1.2*I**0.5)
    # f2_phi = np.exp(-2*I**0.5)

    B0_AA = 0.04166 - 9.6975*1e-5*T
    B0_AC = B0_AA 
    B0_CA = -B0_AA*1.2426
    B1_AC = 0.66838 + 6.081*1e-4*T
    B1_CA = -B1_AC*0.18488
    B2_AC = 7.9545*1e-4
    B2_CA = -B2_AC

    m_A = m_NH3 
    m_C = m_CO2 

    gamma_CO2     = np.exp(2*m_A*(B0_AC) + 3*(m_A**2*B2_AC + 2*m_A*m_C*B2_CA))
    gamma_NH3     = np.exp(2*(m_A*B0_AA + m_C*(B0_CA)) + 3*(2*m_A*m_C*B2_AC + m_C**2*B2_CA))
    
    a_H2O         = np.exp(-M_w*(m_A**2*B0_AA + m_A*m_C*(B0_CA + B0_AC + (B1_AC + B1_CA)) + 6*m_A*m_C*(m_A*B2_AC + m_C*B2_CA)) - np.log(1 + M_w*sum(m)))
    gamma_H2O     = a_H2O/x_H2O

    gamma = [gamma_H2O, gamma_CO2 ,gamma_NH3]

    return gamma

# %%                                                    ############ Gamma CH4 ###########
def fct_gamma_CH4(T,x_CH4):

    a0 = (0,1.360608,0,0.033630,0,0.656974,0,1.763990,0,5.337858,0,-0.024750,0,48.353808,0,-11.580192,0,-0.087295,0,0.558793,0,-23.753020,0,-10.128675,0,-41.212178,0,-31.279868,0,-23.855418,0,-35.675110,0,-33.675110,0,-27.027285,0,-19.026786,0,-37.872252)
    a1 = (0,3.796962,0,-0.703216,0,-12.441339,0,-21.119318,0,-33.298760,0,12.387276,0,17.261174,0,16.384626,0,13.171333,0,13.556732,0,16.573197,0,13.591099,0,5.060082,0,31.289978,0,31.720767,0,37.064849,0,41.544360,0,57.609882,0,54.961702,0,57.204781)
    Tc_CH4 = 190.6 # [K]
    T_R = T/Tc_CH4

    tosum = 0
    for i in range(1,40,2):
        tosum += (a0[i] + a1[i]/T_R)*x_CH4**(0.05+(i-1)/40)
    
    gamma_CH4 = np.exp(tosum)

    return gamma_CH4

# %%
# n_list = [54.815,0.109,0.014,0.681,0.799,0.017,0.092]
# ntot= sum(n_list)
# x = [n/ntot for n in n_list]
# x_H2O, x_CO2, x_NH3, x_HCO3m, x_NH4p, x_CO32m, x_NH2COOm = x

# gamma = gamma_UNIQUAC(333,x_H2O, x_CO2, x_NH3, x_HCO3m, x_NH4p, x_CO32m, x_NH2COOm,species,psi)
# print(gamma)

# gamma_mean = np.exp(np.log(x_H2O) + 1/(ntot)*(sum([n_list[i]*np.log(gamma[i]) for i in range(1,7)])))
# print(gamma_mean)
# print(gamma[0]*x_H2O)
# %%                                                ############## Gamma function V2 ###########

 
def gamma_UNIQUAC_Vext(T, x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm, species, u0_coeff, ut_coeff):
    
# INPUTS : 
# T [K]      : Temperature of liquid phase
# x_i        : Liquid molar fractions
# component  : i-chemical component of interest
# species    : Table containing all the chemical species to consider. !Shall be written with the correct name!
# psi        : psi constants computed from u0 and ut coefficient according to UNIQUAC model

# OUTPUTS : 
# gamma_UNIQUAC : Coefficient of activity of any solute or ionic species computed with UNIQUAC model

# DESCRIPTION : 
# This function computes the coefficient of activity for water, using the extended UNIQUAC model presented
# in Marounina's thesis p34-36. 

                                #####################################################
    
# Extraction of constants :   
    x           = np.array([x_H2O,x_CO2,x_NH3,x_HCO3m,x_Hp,x_NH4p,x_CO32m,x_OHm,x_NH2COOm])

    cmps = range(len(x))

    N           = len(x)
    z           = 10        # Coordination number. Fixed as 10 according to Maurer and Prausnitz (1978)
    r_coeff     = np.array([species['r_coeff'][i] for i in cmps])
    q_coeff     = np.array([species['q_coeff'][i] for i in cmps])
    index_water = getIndexes_onecolumn(species,'H2O')[0]

    N   = 9
    psi = np.zeros((N,N))
    u   = np.zeros((N,N))
        
    for i in range(0,N):
        for j in range(0,N):
            u[i,j]   = u0_coeff.iat[i,j] + ut_coeff.iat[i,j]*(T-298.15)   # u computed with data from litterature (Darde et al. 2010)
        
    for j in range(N):
        for i in range(N):      
            psi[j,i] = np.exp(-(u[j,i]-u[i,i])/T)                         # Computation of psi from UNIQUAC model
    taus = psi

# Computation of necessary constants :

    # xi_ri = np.zeros(N)
    # sum_xi_ri = 0
    # for i in range(N):
    #     xi_ri[i] = x[i]*r_coeff[i]
    #     sum_xi_ri += xi_ri[i]

    # phi = np.array([xi_ri[i]/sum_xi_ri for i in range(N)])
    rsxs   = np.array([r_coeff[i]*x[i] for i in cmps])
    rsxs_sum_inv = 1.0/np.sum(rsxs)
    phis   = np.array([rsxs[i]*rsxs_sum_inv for i in cmps])

    # xi_qi = np.zeros(N)
    # sum_xi_qi = 0
    # for i in range(N):
    #     xi_qi[i] = x[i]*q_coeff[i]
    #     sum_xi_qi += xi_qi[i]
    
    # theta = np.array([xi_qi[i]/sum_xi_qi for i in range(N)])
    qsxs   = np.array([q_coeff[i]*x[i] for i in cmps])
    qsxs_sum_inv = 1.0/np.sum(qsxs)
    thetas = np.array([qsxs[i]*qsxs_sum_inv for i in cmps])

    Ss     = np.array([np.sum([thetas[j]*taus[j][i] for j in cmps]) for i in cmps])
    VsSs   = np.array([thetas[j]/Ss[j] for j in cmps])

# Computation of symmetric coefficients of activity : 

    loggammac = np.zeros(N)
    log_gamma_c = np.zeros(N)

    # for i in range(N):
    #     log_gamma_c[i] = np.log(phi[i]/x[i]) + 1.0 - phi[i]/x[i] - (z/2)*q_coeff[i]*(np.log(theta[i]/phi[i]) + 1.0 - theta[i]/phi[i])

    # log_gamma_r = np.zeros(N)
    # for i in range(N):
    #     sum_theta_psi = 0
    #     for j in range(N):
    #         sum_theta_psi += theta[j]*psi[j][i]
    
    # sum_theta_psi_sum_theta_psi = np.zeros(N)
    # for i in range(N):
    #     for j in range(N): 
    #         sum_theta_psi2 = 0

    #         for k in range(N):
    #             sum_theta_psi2 += theta[k]*psi[k][j]

    #         sum_theta_psi_sum_theta_psi[i] += theta[j]*psi[i][j]/sum_theta_psi2

    #     log_gamma_r[i] = q_coeff[i]*(1.0 - np.log(sum_theta_psi) - sum_theta_psi_sum_theta_psi[i])
    loggammar = np.zeros(N)

    for i in cmps:
        x1 = phis[i]/x[i]
        x2 = phis[i]/thetas[i]
        loggammac[i] = (np.log(x1) + 1.0 - x1 - (z/2)*q_coeff[i]*(np.log(x2) + 1.0 - x2))
        loggammar[i] = (q_coeff[i]*(1.0 - np.log(Ss[i]) - np.sum([taus[i][j]*VsSs[j] for j in cmps])))
    
    q_w  = q_coeff[index_water]
    r_w  = r_coeff[index_water]  
    
# Computation of infinite dilution terms 

    loggammacinf = np.zeros(N)
    loggammarinf = np.zeros(N)

    for i in cmps:
        loggammacinf[i] = (np.log(r_coeff[i]/r_w) + 1 - (r_coeff[i]/r_w) - 0.5*z*q_coeff[i]*(np.log(r_coeff[i]*q_w/(r_w*q_coeff[i])) + 1 - r_coeff[i]*q_w/(r_w*q_coeff[i])))
        loggammarinf[i] = (q_coeff[i]*(1 - np.log(taus[index_water,i]) - taus[i,index_water]))

    # log_gamma_cinf = np.zeros(N)
    # log_gamma_rinf = np.zeros(N)
    # for i in range(N):
    #     log_gamma_cinf[i] = np.log(r_coeff[i]/r_w) + 1 - (r_coeff[i]/r_w) - 0.5*z*q_coeff[i]*(np.log((r_coeff[i]*q_w)/(r_w*q_coeff[i])) + 1 - (r_coeff[i]*q_w)/(r_w*q_coeff[i]))
    #     log_gamma_rinf[i] = q_coeff[i]*(1 - np.log(psi[index_water,i]) - psi[i,index_water])

# Computation of Debye_Huckel term 

    A     = 1.31 + 1.335*1e-3*(T - 273.15) + 1.164*1e-5*((T - 273.15)**2)  # VALIDE BETWEEN [273.15-373.15] K
    Mw    = 18*1e-3                 # Molar masse of water
    b     = 1.5                     # Constant [kg**1/2 mol**-1/2]
    z_dh  = species['charge']
    
    m = np.zeros(N)
    I = 0

    for i in range(N):
        m[i] = x[i]/(Mw*x_H2O)
        I += 0.5*m[i]*z_dh[i]**2    # Ionic strength

    loggammaDH = []
    for i in cmps:
        loggammaDH.append(-(z_dh[i]**2)*(A*I**0.5)/(1 + b*I**0.5))

    sigma = lambda x: 3*(1 + x - 1/(1 + x) - 2*np.log(1 + x))
    loggammaDH_w = 2/3*Mw*A*I**(3/2)*sigma(b*I**(0.5))
    
# Final gamma computation

    # print('c:',np.exp(log_gamma_c), np.exp(loggammac))
    # print('r:',np.exp(log_gamma_r), np.exp(loggammar))
    # print('cinf:',np.exp(log_gamma_cinf), np.exp(loggammacinf))
    # print('rinf:',np.exp(log_gamma_rinf), np.exp(loggammarinf))

    gamma_UNIQUAC = np.zeros(N)
    for i in range(N):
        if i == 0:
            gamma_UNIQUAC[i] = (np.exp(loggammac[i] - loggammacinf[i] + loggammar[i] - loggammarinf[i] + loggammaDH_w))
        else:
            gamma_UNIQUAC[i] = (np.exp(loggammac[i] - loggammacinf[i] + loggammar[i] - loggammarinf[i] + loggammaDH[i])) 

    return gamma_UNIQUAC      

# %%                                               ############## Gamma function V3 ###########

def gamma_UNIQUAC_V3(T,x,species,u0_coeff,ut_coeff):

# INPUTS :
    N = len(x)
    ri = np.array([species['r_coeff'][i] for i in range(N)])
    qi = np.array([species['q_coeff'][i] for i in range(N)])

    tau = np.zeros((N,N))
    u   = np.zeros((N,N))
        
    for i in range(0,N):
        for j in range(0,N):
            u[i,j]   = u0_coeff.iat[i,j] + ut_coeff.iat[i,j]*(T-298.15)   # u computed with data from litterature (Darde et al. 2010)
        
    for j in range(N):
        for i in range(N):      
            tau[j,i] = np.exp(-(u[j,i]-u[i,i])/T)                         # Computation of psi from UNIQUAC model

    rx = np.dot(x, ri)
    qx = np.dot(x, qi)
    phi_x = ri/rx
    tetha = x*qi / qx
    phi_tetha = (ri*qx) / (qi*rx)

    lngamac = np.log(phi_x) + 1. - phi_x
    lngamac -= 5.*qi*(np.log(phi_tetha) + 1. - phi_tetha)

    # residual part
    a0 = u0_coeff.values
    a1 = ut_coeff.values
    Aij = a0 + a1 * T
    tau = np.exp(-Aij/T)
    Sj = np.matmul(tetha, tau)
    SumA = np.matmul(tau, (tetha/Sj))
    lngamar = 1. - np.log(Sj) - SumA
    lngamar *= qi

    # Coefficients at infinite dilution
    # lngamacinf = np.log(ri/ri[0]) + 1. - ri/ri[0] - 5.*qi*(np.log(qi/qi[0]) + 1. - qi/qi[0])
    # lngamarinf = qi*(1. - np.log(tau[0,:]) - tau[:,0])

    #lngama = lngamac + lngamar - lngamacinf - lngamarinf
    lngama = lngamac + lngamar

    return np.exp(lngama)
# %%
# # test
# #
M_H2O = 18.01528*1e-3
m = np.array([54.551,0.788,4.43e-4,0.957,1.01e-6,0.959,7e-9,2.59e-4,7.55e-4])
m_H2O, m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm = m
x_H2O = (1000/M_H2O)/(1000/M_H2O + np.sum(m))
x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm = np.array([m_CO2, m_NH3, m_HCO3m, m_Hp, m_NH4p, m_CO32m, m_OHm, m_NH2COOm])*(x_H2O*M_H2O)
x = np.array([x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm])
#print(gamma_UNIQUAC_Vext(T,x_H2O, x_CO2, x_NH3, x_HCO3m, x_Hp, x_NH4p, x_CO32m, x_OHm, x_NH2COOm,species,u0_coeff,ut_coeff))       
#print(gamma_UNIQUAC_V3(300,x,species,u0_coeff,ut_coeff))
# %%
