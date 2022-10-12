import numpy as np

def mass_to_diameter(mm): # in g
    m = mm/1e3 # to kg
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return 2*((3/(4*np.pi))*m/rho)**(1/3) # km

def diameter_to_mass(D): #km
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return rho*(4/3)*np.pi*(D/2)**3*1e3 # mass in g