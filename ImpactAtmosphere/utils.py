import numpy as np
from . import constants as const

def sat_pressure_H2O(T):
    return 1.0e6*np.exp(const.L_H2O*const.mu_H2O/const.R*(1/373.0 - 1/T)) # dynes/cm^2

# OLR for a steam atmosphere
def OLR(T):
    return 1.6e5 + 500.0*np.maximum(T-1200.0,0.0) # ergs/cm2/s

def mass_to_diameter(mm): # in g
    m = mm/1e3 # to kg
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return 2*((3/(4*np.pi))*m/rho)**(1/3) # km

def diameter_to_mass(D): #km
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return rho*(4/3)*np.pi*(D/2)**3*1e3 # mass in g