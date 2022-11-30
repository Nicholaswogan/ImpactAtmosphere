import numpy as np
from . import constants as const
import numba as nb
from scipy import optimize

def sat_pressure_H2O(T):
    return 1.0e6*np.exp(const.L_H2O*const.mu_H2O/const.R*(1/373.0 - 1/T)) # dynes/cm^2

# net outgoing flux for a steam atmosphere at 4.0 billion years ago
# on Earth. Computed with clima: 
# https://github.com/Nicholaswogan/clima/tree/72bf1781d0cf400e2161c68302ed1b827ee72c89
def net_outgoing_flux(T):
    return 83e3 + 1000.0*np.maximum(T-1750.0,0.0) # ergs/cm2/s

def mass_to_diameter(mm): # in g
    m = mm/1e3 # to kg
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return 2*((3/(4*np.pi))*m/rho)**(1/3) # km

def diameter_to_mass(D): #km
    rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
    return rho*(4/3)*np.pi*(D/2)**3*1e3 # mass in g

@nb.njit
def kc_terms(T, PP, mol_frac):
    """This returns some of the terms in the Kress and Carmichael (1991) relation
    between Fe2O3, FeO and oxygen fugacity

    Parameters
    ----------
    T : float
        Temperature of magma (K)
    PP : float
        Pressure (bars)
    mol_frac : ndarray[double, ndim=1]
        array of mole fractions of species in the magma

    Returns
    -------
    float
        The Kress and Carmichael relation looks like
        `ln(X_Fe2O3/X_FeO) = a*ln(fO2) + tmp`
        This function returns `tmp`
    """    
    P = PP*1e5 # convert from bar to Pa
    
    # constants
    T_0 = 1673
    a, b, c, e, f, g, h = 0.196, 1.1492e4, -6.675, -3.36, -7.01e-7, -1.54e-10, 3.85e-17
    # d = {'Al2O3': -2.243, 'CaO': 3.201, 'Fe2O3': -1.828, 'FeO': -1.828,'K2O': 6.215, 'Na2O': 5.854}
    
    B = b/T
    C = c
    D = (-1.828*mol_frac[0]) + \
        (-1.828*mol_frac[1]) + \
        (-2.243*mol_frac[2]) + \
        (3.201*mol_frac[3]) + \
        (6.215*mol_frac[4]) + \
        (5.854*mol_frac[5])
    E = e*(1 - T_0/T - np.log(T/T_0))
    F = f*P/T
    G = g*(T-T_0)*P/T
    H = h*P**2/T
    
    tmp = B + C + D + E + F + G + H
    return tmp

@nb.njit
def kc_fO2(T, P, mol_frac):
    """The Kress and Carmichael (1991) relation, solved for oxygen fugacity

    Parameters
    ----------
    T : float
        Temperature of magma (K)
    PP : float
        Pressure (bars)
    mol_frac : ndarray[double, ndim=1]
        array of mole fractions of species in the magma

    Returns
    -------
    float
        Oxygen fugacity in bars
    """    
    a = 0.196
    tmp = kc_terms(T, P, mol_frac)
    # Fe2O3 is index 0, and FeO is index 1
    fO2 = np.exp((np.log(mol_frac[0]/mol_frac[1]) - tmp)/a)
    return fO2

@nb.njit
def alter_melt_composition_objective(x, mols_init, T, P, fO2):
    """This is the objective function for `alter_melt_composition_to_fO2`
    """
    Fe2O3, FeO = 10.0**x # mols/g
    mols = mols_init.copy()
    mols[0] = Fe2O3
    mols[1] = FeO
    
    Fe_tot = (0.5*mols_init[0] + mols_init[1])
    Fe_tot_new = (0.5*Fe2O3 + FeO)
    
    mol_frac = mols/np.sum(mols)
    fO2_eq = kc_fO2(T, P, mol_frac)
    
    res = np.array([
        np.log10(fO2_eq) - np.log10(fO2),
        Fe_tot - Fe_tot_new
    ])
    return res

@nb.njit
def log10fO2_FMQ(T):
    """Oxygen fugacity at FMQ from Zahnle et al. 2020"""
    f_O2_QFM = 3.015e-4*T**3.449*np.exp(-53649.0/T) * 1.01325 # conversion from atm to bars
    return np.log10(f_O2_QFM) # log10(bars)

def alter_melt_composition_to_fO2(mols, T, P, DFMQ):
    """Alters a melt composition so that it has a certain redox state."""
    log10_FMQ = log10fO2_FMQ(T)
    fO2 = 10.0**(log10_FMQ + DFMQ) # target fO2
    
    init = np.log10([mols[0], mols[1]])
    sol = optimize.root(alter_melt_composition_objective, init, args = (mols, T, P, fO2))
    if not sol.success:
        raise Exception("root solve failed")

    mols_eq = mols.copy()
    mols_eq[0] = 10.0**sol.x[0]
    mols_eq[1] = 10.0**sol.x[1]
    
    return mols_eq

def alter_basalt_melt_to_fO2(T, P, DFMQ):
    basalt = BasaltComposition()
    return alter_melt_composition_to_fO2(basalt.mols, T, P, DFMQ)

@nb.njit
def rates_melt_reactions(T, P, n_Fe2O3, n_FeO, fO2, mols_melt):
    """Computes the rates of change of O2, FeO and Fe2O3 from the 
    reactions 0.5O2 + 2 FeO <=> Fe2O3. The rates
    are made up to be fast and to push the reaction toward an equilibrium
    state, such that it satisfies the Kress and Carmichael (1991) relation.

    Parameters
    ----------
    T : float
        Temperature in K
    P : float
        Pressure in bar
    n_Fe2O3 : float
        mol Fe2O3 / g melt.
    n_FeO : float
        mol FeO / g melt
    mols_melt : ndarray[double, ndim=1]
        mol/g in the melt of interest. Order of species is given in
        BasaltComposition class.

    Returns
    -------
    tuple
        Rates of change of O2, FeO and Fe2O3 in mol/(g*s).
    """    
    mols = mols_melt.copy()
    mols[0] = n_Fe2O3 # mols/g
    mols[1] = n_FeO # mols/g
    mol_frac = mols/np.sum(mols)

    # forward reaction rate. Totally made up. Chosen to be fast.
    k_forward = 1e-10
    # This is an equilibrium constant for the reaction
    # as determined by Kress and Carmichael (1991)
    Keq = np.exp(kc_terms(T, P, mol_frac))
    # reverse rate constant, so that rates push the reaction
    # toward equilbrium
    k_reverse = k_forward/Keq 

    a = 0.196 # K+C constant
    dn_O2_dt = 0.5*(k_reverse*n_Fe2O3 - k_forward*n_FeO*fO2**a) # mol O2/(g melt*s)
    dn_FeO_dt = 2*(k_reverse*n_Fe2O3 - k_forward*n_FeO*fO2**a) # mol FeO/(g melt*s)
    dn_Fe2O3_dt = -k_reverse*n_Fe2O3 + k_forward*n_FeO*fO2**a # mol Fe2O3/(g melt*s)

    return dn_O2_dt, dn_FeO_dt, dn_Fe2O3_dt

@nb.njit
def rates_melt_H2O_reaction(P_H2O, n_H2O_melt):
    """Computes rates of the reaction
    H2O(g) <=> H2O(aq)
    Where H2O(g) is gas-phase and H2O(aq) is H2O dissolved in magma.
    Here I use a solubility expression from Itcovitz et al. (2022) (Eq. 19)

    Parameters
    ----------
    P_H2O : float
        Partial pressure of H2O in bar
    n_H2O_melt : float
        concentration of H2O in melt (mol / g melt)

    Returns
    -------
    tuple
        Rate of change of H2O in the gas and melt in mol H2O / (g melt * s)
    """        
    # conver to mass fraction (g H2O / g melt)
    X_H2O_melt = n_H2O_melt*const.mu_H2O

    # forward rate constant (made up)
    kf = 1.0e-10
    Keq = 0.00021503488089144967 # = 6.8e-8*(10^5)^0.7
    kr = kf/Keq

    dn_H2O_gas_dt = -kf*P_H2O**0.7 + kr*X_H2O_melt
    dn_H2O_melt_dt = kf*P_H2O**0.7 - kr*X_H2O_melt

    return dn_H2O_gas_dt, dn_H2O_melt_dt

# Composition of basalt
class BasaltComposition():
    def __init__(self):
        # species in the basalt (we ignore H2O, and assume it is minor)
        self.species = [
            'Fe2O3', 
            'FeO', 
            'Al2O3', 
            'CaO', 
            'K2O', 
            'Na2O',
            'SiO2',
            'MnO',
            'MgO',
            'TiO2',
            'P2O5'
        ]
        # mass fractions (g/g)
        self.mass_frac = np.array([
            0.0318, # Fe2O3
            0.0685, # FeO
            0.1318, # Al2O3
            0.1234, # CaO
            0.0093, # K2O
            0.0218, # Na2O
            0.4995, # SiO2
            0.0015, # MnO
            0.0998, # MgO
            0.0101, # TiO2
            0.0025 # P2O5
        ])
        # molar masses (g/mol)
        self.mu = np.array([
            160.,  
            72., 
            102.,  
            56.,  
            94.,  
            62.,  
            60.,  
            71.,  
            40.,  
            80., 
            142.
        ])

        # mols/(g melt)
        self.mols = np.empty(self.mass_frac.shape[0])
        for i in range(self.mols.shape[0]):
            self.mols[i] = self.mass_frac[i]/self.mu[i]

        # mol fraction
        self.mol_frac = self.mols/np.sum(self.mols)