import numpy as np
from scipy.integrate import solve_ivp
import sys
from EvolveAtmFort import atmos
atmos.setup()
from EvolveAtmFort import diffusion


def integrate(tspan,Ninit_dict,out_dict=True,**kwargs):
    '''
    Evolves a Hadean Earth atmosphere using a simple 0-D photochemical
    model from time tspan[0] to tspan[1] given the initial conditions Ninit_dict.
    Uses `scipy.integrate.solve_ivp`.

    Parameters
    ----------
    tspan : list
        Will integrate from time tspan[0] to tspan[1] in seconds.
    Ninit_dict : dict
        Dictionary containing initial atmosphere in molecules/cm^2. Must contain
        the following 9 dictionary items: H2, CO, CO2, CH4, N2, NH3, HCN, C2Hn, Haze
    out_dict: bool
        If true, then the output will be a dictionary. If false, then the output
        will be a numpy array containing the abundance of each molecule at each
        timestep in molecules/cm^2.
    **kwargs:
        The same optional arguments as scipy.integrate.solve_ivp,
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

    Returns
    -------
    out: dict
        Dictionary defining the state of the atmosphere as a function of time.
    '''

    # species
    species = ['H2O','H2','CO','CO2','CH4','N2','NH3', 'HCN', 'C2Hn', 'Haze']

    # dict to array
    Ninit = np.zeros(len(species))
    Ninit[0] = 1. # H2O doesn't matter
    Ninit[1:] = np.array([Ninit_dict[key] for key in species[1:]])

    # set uv to 100. then integrate
    atmos.tau_uv_init = 100.
    sol = solve_ivp(atmos.rhs,tspan,Ninit,**kwargs)

    # only take positive solution (it should be almost positive)
    N_vals = np.clip(sol.y.T,0,np.inf)
    tvals = sol.t

    # check if solver succeeded
    if not sol.success:
        sys.exit('ODE integration failed')

    if not out_dict:
        return N_vals
    elif out_dict:
        # re-run the RHS of the ODEs but this time save interesting stuff
        dNdt = np.zeros(N_vals.shape)
        pressure = np.zeros(len(tvals))
        T_surf = np.zeros(len(tvals))
        tau_uv = np.zeros(len(tvals))
        N_H2O = np.zeros(len(tvals))
        mubar = np.zeros(len(tvals))
        atmos.tau_uv_init = 100.
        for i in range(0,len(tvals)):
            dNdt[i],pressure[i],T_surf[i],tau_uv[i],N_H2O[i],mubar[i] = \
            atmos.rhs_verbose(tvals[i],N_vals[i])

        # build a dictionary of output
        out = {}
        out['Ntot'] = np.sum(N_vals[:,:-3],axis=1)-\
                      N_vals[:,species.index('H2O')] + N_H2O
        out['press'] = pressure
        out['Tsurf'] = T_surf
        out['tau_uv'] = tau_uv
        out['mubar'] = mubar
        out['time'] = tvals
        out['dNHCN_dt'] = dNdt[:,species.index('HCN')]
        out['dNHaze_dt'] = dNdt[:,species.index('Haze')]
        out['H2O'] = N_H2O/out['Ntot']
        out['NH2O_strat'] = N_vals[:,species.index('H2O')]
        for i in range(len(species)):
            if species[i] != 'H2O':
                out[species[i]] = N_vals[:,i]/out['Ntot']

        return out

def diffuse(PhiHCN, Ts = 298, Ps = 1, mubar = 28.0, pH = 7, \
            Kzz = 1.0e5, top_atm = 60.0e5, nz = 60, T_trop = 180, \
            P_trop = 0.1, **kwargs):
    '''
    Calculates the HCN mixing ratio as a function of altitude for a
    given HCN production rate (PhiHCN)

    Parameters
    ----------
    PhiHCN: float
        The HCN production rate (molecules/cm^2/s)
    Ts: float, optional
        The surface temperature (K)
    Ps: float, optional
        The surface pressure (bar)
    mubar: float, optional
        The mean molar weight of the atmosphere (g/mol)
    pH: float, optional
        The pH of the ocean (cm/s)
    Kzz: float, optional
        The eddy diffusion coefficient (cm^2/s)
    top_atm: float, optional
        The top of the atmosphere (cm)
    nz: integer, optional
        The number of vertical descritization in the atmosphere.
    T_trop: float, optional
        Tropospheric temperature (K)
    P_trop: float, optional
        Tropospheric pressure (bar)

    Returns
    -------
    alt: numpy array length nz
        The altitude of each HCN mixing ratio (cm).
    fHCN: numpy array length nz
        The HCN mixing ratio as a function of altitude in the atmosphere.
    '''

    vd_HCN = HCN_vdep(Ts,pH,**kwargs)

    alt, fHCN = diffusion.diffuse(PhiHCN, Ts, Ps, mubar, vd_HCN, \
                                  Kzz, top_atm, nz, T_trop, P_trop)

    return alt, fHCN

def HCN_vdep(T, pH, vo = 1.2e-5, zs = 100e2, zd = 4000e2):
    '''
    Calculates deposition velocity (cm/s) of HCN into the ocean assuming
    it is destroyed from hydrolysis.

    Parameters
    ----------
    T: float
        Temperature of the ocean (L)
    pH: float
        The pH of the ocean (unit-less)
    vo: float, optional
        The turnover velocity of the ocean (cm/s)
    zs: float, optional
        The depth of the surface ocean (cm)
    zd: float, optional
        The depth of the deep ocean (cm). Note zs + zd = total ocean depth.

    Returns
    -------
    vd_HCN: float
        Deposition velocity of HCN into the ocean (cm/s)
    '''
    vp_HCN = 0.005 # piston velocity of HCN (cm/s). WARNING! need to confirm this is reasonable
    CC = 6.022e20 # avagadros number/1e3 (molecules/mol)
    k = 1.3807e-16 # boltzman (cgs units)

    # HCN henry's law constant
    lnkH=8205.7/T-25.323   #in M/atm
    alpha_HCN=np.exp(lnkH)/1.013 # conver to M/bar

    ktot = HCN_hydrolysis_rate(T,pH)

    vd_HCN = k*T*1e-6*alpha_HCN*CC*ktot*vp_HCN*(ktot*zd*zs+vo*(zd+zs))\
                    /(ktot*zd*(vp_HCN+ktot*zs)+vo*(vp_HCN+ktot*(zd+zs)))
    return vd_HCN

def HCN_hydrolysis_rate(T, pH):
    '''
    Calculates HCN hydrolysis rate following Miyakawa et al. 2001.

    Parameters
    ----------
    T: float
        Temperature of the water (K)
    pH; float
        The pH of the water (unit-less)

    Returns
    -------
    ktot: float
        The HCN hydrolysis rate (1/s)
    '''
    H=10.**(-1.*pH)
    OH=1.0e-14/H

    # acid catalyzed: (in M^-1 s^-1)
    logk1H=-4950./T + 8.43
    k1H=10.**(logk1H)
    # base catalyzed: (in M^-1 s^-1)
    logk1OH=-4240./T + 11.1
    k1OH=10.**(logk1OH)

    pKw=-6.0846+4471.33/T+.017053*T  #from Stribling and Miller 1987, from Miyakawa's original sources (Robinson and Stokes 1959 and Schlesinger and Miller 1973, resp.)
    pKa_HCN=-8.85+3802./T+.01786*T
    Kw=10.**(-1.*pKw)
    Ka_HCN=10.**(-1.*pKa_HCN)

    ktot=k1H*H+(k1OH*Kw/(H+Ka_HCN))  #in s^-1

    return ktot # in s^-1
