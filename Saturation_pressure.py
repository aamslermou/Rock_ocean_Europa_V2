import numpy as np
# This function computes the water saturation pressure in Pa
# VALID BETWEEN 255.9 - 373 K (Antoine law) Stull, 1947 
# VALID BETWEEN 379-573 K : Liu and Lindsay, 1970
# source :(https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=4&Type=ANTOINE&Plot=on#ANTOINE)

def P_sat_water(T):

    if T>373:
        A = 3.55959
        B = 643.748
        C = -198.043
    else:
        A = 4.6543
        B = 1435.264
        C = - 64.848

    PsatH2O = 10**(A - (B/(T+C)))       #[bar]

    return PsatH2O*1e5

# This function the carbon dioxide saturation pressure in Pa
# Valid between 154.15 - 304.15 K
# source : https://www.cheric.org/research/kdb/hcprop/showcoef.php?cmpid=1943&prop=PVP

def P_sat_CO2(T):

    A = -2.403761E+01
    B = -7.062404E+03
    C = 1.663861E+02
    D = 3.368548E-05

    PsatCO2 = np.exp(A*np.log(T) + B/T + C + D*T**2) # [kPa]

    return PsatCO2*1e3

