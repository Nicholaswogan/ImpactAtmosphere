import numpy as np
from . import constants as const

def sat_pressure_H2O(T):
    return 1.0e6*np.exp(const.L_H2O*const.mu_H2O/const.R*(1/373.0 - 1/T)) # dynes/cm^2

# net outgoing flux for a steam atmosphere at 4.0 billion years ago
# on Earth. Computed with clima: 
# https://github.com/Nicholaswogan/clima/tree/72bf1781d0cf400e2161c68302ed1b827ee72c89
def net_outgoing_flux(T):
    return 83e3 + 1000.0*np.maximum(T-1750.0,0.0) # ergs/cm2/s

def oxygen_fugacity_QFM(T, DQFM):
    """From Zahnle et al. 2020"""
    f_O2_QFM = 3.015e-4*T**3.449*np.exp(-53649.0/T) * 1.01325 # conversion from atm to bar
    # in bar
    log10_FO2 = np.log10(f_O2_QFM)
    f_O2 = 10.0**(log10_FO2 + DQFM) 
    return f_O2 # bars

def mass_to_diameter(mm): # in g
    m = mm/1e3 # to kg
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return 2*((3/(4*np.pi))*m/rho)**(1/3) # km

def diameter_to_mass(D): #km
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return rho*(4/3)*np.pi*(D/2)**3*1e3 # mass in g