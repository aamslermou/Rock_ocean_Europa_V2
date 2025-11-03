#%%
import numpy as np
import pandas as pd
#%%                                    #################################### DATA ################################################

                                    
species = pd.DataFrame(np.array([['H2O', 647.3, 22.0, 0.3434, 18.01528, 0, 0.92, 1.4000],
                                 ['CO2', 304.19, 7.383, 0.224,  44.0095, 0, 0.75, 2.45],
                                 ['NH3', 405.6, 11.28, 0.25, 17.03052, 0, 1.6292, 2.9852]]),
                       
                        columns = ['formula', 'Tc', 'Pc', 'w', 'M', 'charge', 'r_coeff', 
                           'q_coeff'])

k_coeff = np.array([[0,  0.1896, 0],
                    [0.1896, 0, 0],
                    [0,0,0]])


#%%                                    ####################### Function to compute the fugacity ##############################

from search_function import getIndexes_onecolumn

def fugacity_coeff_vapor_PG(T,P,component,species,y_H2O,y_CO2,y_NH3,b,A,A_mix,B_mix):

# INPUTS : 
# T [K]     : Temperature of vapor phase
# P [MPa]   : Pressure of the vapor phase
# component : i-chemical component of interest
# species   : Table containing all the chemical species to consider. !Shall be written with the correct name!
# y         : Gaseous molar fractions y_i
# b         : b coefficient from Peng-Robinson equation of state
# A         : A coefficient from Peng-Robinson equation of state
# A_mix     : Mix law coefficient computed from Peng-Robinson equation of state coefficients
# B_mix     : Mix law coefficient computed from Peng-Robinson equation of state coefficients

# OUTPUTS : 
# phi_vap   : Fugacity coefficient 

# DESCRIPTION : 
# This function computes the fucacity coefficient for ONE component in VAPOR phase. 
# Based on the assumption the the vapor phase is a simple mix. We use the Peng-Robinson equation of state
# and the expression on phi_i form Marounina thesis eq 2.23, p31

                                #####################################################
    
# Extraction of constants :
    species = species[0:3]          # To only consider H2O, CO2, NH3
    R       = 8.314
    index   = getIndexes_onecolumn(species,component)[0]
    r1      = -1 - 2**0.5
    r2      = -1 + 2**0.5
    y       = [y_H2O, y_CO2, y_NH3]
    k_index = k_coeff[component]

    a_sum = np.zeros([3,3])

    for i in range(0,3):
        for j in range(0,3):
            a_sum[i] = y[i] = (1-k_coeff.iat[i+1,j+1])*(A[i]*A[j])**0.5 

    a_sum_sum = np.zeros(3)
    for i in range(0,3):
        a_sum_sum[i] = sum(a_sum[:,i])
         
    delta = 2*a_sum_sum[index]/A_mix
# Computation of compressibility factor Z, solving Peng-Robinson equation of state :

    A  = (A_mix*P)/(R**2*T**2)    # Using polynominal form
    B  = (B_mix*P)/(R*T)
    #print(A)
    #print(B)

    m0 = - (A *B - B**2 - B**3)
    m1 = A - 3 * B**2 - 2 * B
    m2 = - (1 - B)

    coeff = [1, m2, m1, m0]
    roots = np.roots(coeff) 
    #print(roots)
    Z     = np.max(roots)
    #print(Z)
# Final computation : 
    bi = b[index]
    #print(v)
    phi_vap = np.exp((bi/B_mix)*(Z-1) - np.log(Z - B) \
                 - A/(B*(r1-r2))*np.log((Z-B*r1)/(Z-B*r2)))*(delta - bi/B_mix) 
    #phi_vap = np.exp((bi/B_mix)*(Z-1) - np.log(P*(v-B_mix)/(R*T)) \
                 #+ A_mix/(R*T*B_mix*(r1-r2))*np.log((v-B_mix*r1)/(v-B_mix*r2))*(delta - bi/B_mix))
    
    return np.real(phi_vap)       

# %%                                   ####################### Function to compute the fugacity v2 ##############################
                                  

def FUG_COEFF_PRG_EOS(R, T, P, species, k_coeff, y): 
  
  N = len(y)
  # components = species  
  Tc = np.array([species['Tc'][i] for i in range(N)])
  Pc = np.array([species['Pc'][i] for i in range(N)])
  w  = np.array([species['w'][i] for i in range(N)])

# a_k, b_k
  a = np.zeros([N])

  b = np.zeros([N])

  for i in range(N):
   
    # Calculate the alpha term
    kappa = 0.37464 + 1.54226 * w[i] - 0.26992 * w[i]**2
    alpha = (1 + kappa * (1 - np.sqrt(T / Tc[i])))**2
    # Calculate the b coefficient using the Peng-Robinson equation of state
    b[i] = 0.0778 *(R * Tc[i]) / Pc[i]  # b coefficient
    # Calculate the a[i] term using the alpha term
    a[i] = alpha * 0.45724 * (R**2 * Tc[i]**2) / Pc[i]
# aij

  A_ij = np.zeros([N,N])

  for i in range(N):    
    for j in range(N):      
      A_ij[i,j] = np.sqrt(a[i]*a[j]) * (1 - k_coeff[i,j])

# A_mix

  A_mix = 0

  xixjaij = np.zeros([N,N])
  for i in range(N):
    for j in range(N):
      xixjaij[i,j] = y[i] * y[j] * A_ij[i,j]
      #A_mix += y[i] * y[j] * A_ij[i,j]
  A_mix = sum(sum(xixjaij))
  sum_xa_Aij = np.zeros([N])


# xi aik

#   xa_ik = np.zeros([N,N])

#   for i in range(N):
    
#     for j in range(N):
#       xa_ik[i,j] = y[i] * A_ij[i,j]

#   xa_ik_sum = np.zeros([N])

#   for i in range(N):
#     xa_ik_sum[i] = sum(xa_ik[:,i])

# # xi xj aij

#   xxa_ij = np.zeros([N,N])

#   for i in range(N):
    
#     for j in range(N):
      
#       xxa_ij[i,j] = y[i] * y[j] * A_ij[i,j]

#   A_mix = np.zeros([N])

#   for i in range(N):
#     A_mix[i] = sum(xxa_ij[:,i])

#   A_mix = sum(A_mix)

# b(T)

  B_mix = np.zeros([N]) 

  for i in np.arange(N):
      B_mix[i] = (b[i] * y[i])

  B_mix = sum(B_mix)
# Peng Robinson coefficients

  A = A_mix * P / (R**2 * T**2)
  B = B_mix * P / (R * T)

  m0 = - (A *B - B**2 - B**3)
  m1 = A - 3 * B**2 - 2 * B
  m2 = - (1 - B)

  coeff = [1, m2, m1, m0]
  Z = np.roots(coeff)    
  Z = np.real(Z)
  Z = np.max(Z)
# Computing phi

  phi_V_list = np.zeros([N])

  for i in range(N):
      for j in range(N):
        sum_xa_Aij[i] += y[j] * A_ij[j,i]
      phi_V_list[i] = np.exp( (b[i]/B_mix)*(Z-1) - np.log(Z - B) - (A / (2 * 2**0.5 * B)) * \
                             ((2*sum_xa_Aij[i]/A_mix) - (b[i]/B_mix)) * np.log((Z + 2.414 * B) / (Z - 0.414 * B)))
      #phi_V = np.exp(b[i] / B_mix *(np.max(Z) - 1) - np.log(np.max(Z) - B) - A / (2 * 2**.5 * B) * \
      #               (2 * xa_ik_sum[i] / A_mix - b[i] / B_mix)  * \
      #                np.log((np.max(Z) + 2.414 * B) / (np.max(Z) - 0.414 * B))) 
      #print(phi_V) 
      #phi_V_list[i] = phi_V
      
#   c_fug_L_tod = np.zeros([N])

#   for i in np.arange(N):
#     c_fug_L = np.exp(b_k[i] / b_T *(np.min(odin) - 1) - np.log(np.min(odin) - B) - A / (2 * 2**.5 * B) * (2 * xa_ik_sum[i] / a_T - b_k[i] / b_T)  * np.log((np.min(odin) + 2.414 * B) / (np.min(odin) - 0.414 * B)))  
#     c_fug_L_tod[i] = c_fug_L
    
  return phi_V_list, np.max(Z)#, c_fug_L_tod

#%%                                                     ########### Test function fugacity ##########

# components = species[0:3]
# y_H2O = 0.1
# y_CO2 = 0.9
# y_NH3 = 0.0

# y0 = [y_H2O,y_CO2,y_NH3]
# T     = 333         # [K]
# P0    = 2           # [MPa]
# R     = 8.314
# # T_c   = components['Tc']
# # P_c   = components['Pc']
# # w     = components['w']

# # a     = (0.457235*R**2*T_c**2)/P_c
# # kappa = 0.37464 + 1.54226*w - 0.26992*w**2
# # alpha = (1 + kappa*(1 - (T/T_c)**0.5))**2
# # A     = a*alpha
# # b     = 0.077796*R*T_c/P_c

# P_list    = np.linspace(2,40,50)

# phi_H2O_list = []
# phi_CO2_list = []
# phi_NH3_list = []

# for P in P_list:
# #     A_mix = 0    
# #     for i in range(0,3):
# #         for j in range(0,3):
# #             A_ij   = (1-k_coeff.iat[i+1,j+1])*(A[i]*A[j])**0.5 
# #             A_mix += y0[i]*y0[j]*A_ij
# # #print(A_mix)

# #     B_mix = 0
# #     for i in range(0,3):
# #         B_mix += y0[i]*b[i]  

#     # Phi_H2O = fugacity_coeff_vapor_PG(T,P,'H2O',components,y_H2O,y_CO2,y_NH3,b,A,A_mix,B_mix)
#     # Phi_CO2 = fugacity_coeff_vapor_PG(T,P,'CO2',components,y_H2O,y_CO2,y_NH3,b,A,A_mix,B_mix)
#     # Phi_NH3 = fugacity_coeff_vapor_PG(T,P,'NH3',components,y_H2O,y_CO2,y_NH3,b,A,A_mix,B_mix)
#     Phi = FUG_COEFF_PRG_EOS(8.314, T, P, species, k_coeff, y0)
#     phi_H2O_list.append(Phi[0])
#     phi_CO2_list.append(Phi[1])
#     phi_NH3_list.append(Phi[2])

# multiple = 10
# P_bar    = [x*multiple for x in P_list]

# plt.figure(1)

# plt.plot(P_bar,phi_H2O_list, color='red', label='$Phi_{H2O}$')
# plt.plot(P_bar,phi_CO2_list, color='green', label='$Phi_{CO2}$')
# plt.plot(P_bar,phi_NH3_list, color='blue', label='$Phi_{NH3}$')

# plt.title('T = 300 K')
# plt.xlabel('Pression [bar]')
# plt.ylabel('Phi')
# plt.legend()

# plt.show()
#%%                                                     ############ Fugacity 2 species ###########

from search_function import getIndexes_onecolumn


def fugacity_coeff_vapor_PG_2species(T,P,component,species,y_H2O,y_CO2,b,A,A_mix,B_mix):

# INPUTS : 
# T [K]     : Temperature of vapor phase
# P [MPa]   : Pressure of the vapor phase
# component : i-chemical component of interest
# species   : Table containing all the chemical species to consider. !Shall be written with the correct name!
# y         : Gaseous molar fractions y_i
# b         : b coefficient from Peng-Robinson equation of state
# A         : A coefficient from Peng-Robinson equation of state
# A_mix     : Mix law coefficient computed from Peng-Robinson equation of state coefficients
# B_mix     : Mix law coefficient computed from Peng-Robinson equation of state coefficients

# OUTPUTS : 
# phi_vap   : Fugacity coefficient 

# DESCRIPTION : 
# This function computes the fucacity coefficient for ONE component in VAPOR phase. 
# Based on the assumption the the vapor phase is a simple mix. We use the Peng-Robinson equation of state
# and the expression on phi_i form Marounina thesis eq 2.23, p31

                                #####################################################
    
# Extraction of constants :
    species = species[0:2]          # To only consider H2O, CO2
    R       = 8.314
    index   = getIndexes_onecolumn(species,component)[0]
    r1      = -1 - 2**0.5
    r2      = -1 + 2**0.5
    y       = [y_H2O, y_CO2]
    k_index = k_coeff[component]

    delta = 2*A[index]/A_mix * np.sum(y*(A**0.5)*(1-k_index))

# Computation of compressibility factor Z, solving Peng-Robinson equation of state :

    A  = (A_mix*P)/((R*T)**2)    # Using polynominal form
    B  = (B_mix*P)/(R*T)
    #print(A)
    #print(B)

    m0 = B**3 + B**2 - A*B
    m1 = - 3*B**2 - 2*B + A
    m2 = B - 1

    coeff = [1, m2, m1, m0]
    roots = np.roots(coeff) 
    #print(roots)
    Z     = np.max(roots)
    #print(Z)
# Final computation : 
    bi = b[index]
    v  = R*T*Z/P
    #print(v)
    
    phi_vap = np.exp((bi/B_mix)*(Z-1) - np.log(P*(v-B_mix)/(R*T)) \
                 + A_mix/(R*T*B_mix*(r1-r2))*np.log((v-B_mix*r1)/(v-B_mix*r2))*(delta - bi/B_mix))

    return np.real(phi_vap)

# %%

# species_f = species[0:2]          # To only consider H2O, CO2
# T_c       = species_f['Tc']
# P_c       = species_f['Pc']
# w         = species_f['w']
# R         = 8.314
# T         = 300 
# P_list    = np.linspace(20*1e5,300*1e5,50)

# phi_H2O_list = []
# phi_2_list   = []

# for P in P_list:
#     a     = (0.457235*R**2*T_c**2)/P_c
#     kappa = 0.37464 + 1.54226*w - 0.26992*w**2
#     alpha = (1 + kappa*(1 - (T/T_c)**0.5))**2
#     A     = a*alpha
#     b     = 0.077796*R*T_c/P_c

#     y0    = [0.1,0.9]

#     A_mix  = 0    
#     for i in range(0,2):
#         for j in range(0,2):
#             A_ij   = (1-k_coeff.iat[i+1,j+1])*(A[i]*A[j])**0.5
#             A_mix += y0[i]*y0[j]*A_ij
 
#     B_mix = 0
#     for i in range(0,2):
#         B_mix += y0[i]*b[i]   
     
#     phi_H2O = fugacity_coeff_vapor_PG_2species(T,P,'H2O',species,y0[0],y0[1],b,A,A_mix,B_mix)
#     phi_2   = fugacity_coeff_vapor_PG_2species(T,P,'CO2',species,y0[0],y0[1],b,A,A_mix,B_mix)    

#     phi_H2O_list.append(phi_H2O)
#     phi_2_list.append(phi_2)

# multiple = 1e-5
# P_bar    = [x*multiple for x in P_list]

# plt.figure(2)

# plt.plot(P_bar,phi_H2O_list, color='red', label='$Phi_{H2O}$')
# plt.plot(P_bar,phi_2_list, color='green', label='$Phi_{CO2}$')

# plt.title('T = 300 K')
# plt.xlabel('Pression [bar]')
# plt.ylabel('Phi')
# plt.legend()

# plt.show()



# %%
#print(FUG_COEFF_PRG_EOS(8.314, 300, 3, species, k_coeff, [0.0007176439498046115, 0.999282342313496, 1.3736699445298989e-08]))

# %%
