import numpy as np
from scipy.integrate import solve_ivp
import sys
from EvolveAtmFort import atmos
atmos.setup()
from EvolveAtmFort import diffusion


class EvolveAtm:

    def __init__(self):
        self.species = ['H2O','H2','CO','CO2','CH4','N2','NH3', 'HCN', 'C2Hn', 'Haze']
        self.atmos = atmos
        self.diffusion = diffusion

    def rhs(self,t,y):
        '''
        Returns the right-hand-side of the ODEs defining
        change in the atmosphere
        '''
        return self.atmos.rhs(t,y)

    def integrate(self,tspan,Ninit_dict,out_dict=True,**kwargs):
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
        '''

        # dict to array
        Ninit = np.zeros(len(self.species))
        Ninit[0] = 1. # H2O doesn't matter
        Ninit[1:] = np.array([Ninit_dict[key] for key in self.species[1:]])

        # set uv to 100. then integrate
        self.atmos.tau_uv_init = 100.
        sol = solve_ivp(self.atmos.rhs,tspan,Ninit,**kwargs)

        # only take positive solution
        N_vals = np.clip(sol.y.T,1,np.inf)
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
            self.atmos.tau_uv_init = 100.
            for i in range(0,len(tvals)):
                dNdt[i],pressure[i],T_surf[i],tau_uv[i],N_H2O[i],mubar[i] = \
                self.atmos.rhs_verbose(tvals[i],N_vals[i])

            # build a dictionary of output
            out = {}
            out['Ntot'] = np.sum(N_vals[:,:-3],axis=1)-\
                          N_vals[:,self.species.index('H2O')] + N_H2O
            out['press'] = pressure
            out['Tsurf'] = T_surf
            out['tau_uv'] = tau_uv
            out['mubar'] = mubar
            out['time'] = tvals
            out['dNHCN_dt'] = dNdt[:,self.species.index('HCN')]
            out['dNHaze_dt'] = dNdt[:,self.species.index('Haze')]
            out['H2O'] = N_H2O/out['Ntot']
            out['NH2O_strat'] = N_vals[:,self.species.index('H2O')]
            for i in range(len(self.species)):
                if self.species[i] != 'H2O':
                    out[self.species[i]] = N_vals[:,i]/out['Ntot']

            return out

    def diffuse(self,PhiHCN, Ts = 280, Ps = 1, mubar = 28.0, vd_HCN = 5.0e-3, \
                Kzz = 1.0e5, top_atm = 60.0e5, nz = 60, T_trop = 180, \
                P_trop = 0.1):
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
        vd_HCN: float, optional
            The deposition velocity of HCN into the ground (cm/s)
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
        alt: numpy array
            The altitude of each HCN mixing ratio (cm).
        fHCN: numpy array
            The HCN mixing ratio as a function of altitude in the atmosphere.
        '''

        alt, fHCN = self.diffusion.diffuse(PhiHCN, Ts, Ps, mubar, vd_HCN, \
                                           Kzz, top_atm, nz, T_trop, P_trop)

        return alt, fHCN

    def HCN_vdep(self,k1,vp_HCN = 0.005, alpha_HCN = 7.5, vo = 2e-5, \
                 zs = 100e2, zd = 4000e2, T = 298):
        '''
        Calculates deposition velocity of HCN into the ocean assuming
        it is destroyed from hydrolysis.
        '''
        CC = 6.02e20
        k = 1.3807e-16 #cgs

        vd_HCN = k*T*1e-6*alpha_HCN*CC*k1*vp_HCN*(k1*zd*zs+vo*(zd+zs))\
                 /(k1*zd*(vp_HCN+k1*zs)+vo*(vp_HCN+k1*(zd+zs)))
        return vd_HCN

    def HCN_hydrolysis_rate(self):
        '''
        Calculates HCN hydrolysis rate.
        '''
        return None
