import numpy as np
from scipy.integrate import solve_ivp
import sys

from EvolveAtmFort import atmos
atmos.setup()
from .hydrolysis_rates import HCN_hydrolysis_rate, HCONH2_hydrolysis_rate
from .EvolveAtm import diffuseHCN
from .SteamAtm import SteamAtm


class LakeAfterImpact:
    """This class contains a 0-D photochemical model coupled to a surface lake.
    """

    def __init__(self,T_surf,pH_ocean,T_lake,pH_lake,Kzz = 1.e5,\
                   HCNalt = 60.0e5, nz = 60, T_trop = 180., \
                   P_trop = 0.1, vo = 1.2e-5, zs = 100.e2, zd = 4000.e2, \
                   fH2O_trop = 1.0e-6):
        """Set parameters used for the simulation.

        Parameters
        ----------
        T_surf : float
            Surface temperature of the atmosphere and temperature
            of the ocean (K).
        pH_ocean : float
            pH of the ocean (no units).
        T_lake : float
            Temperature of the small lake on land (K).
        pH_lake : float
            pH of the small lake (unit-less)
        Kzz : float, optional
            Eddy diffusion coefficient (cm2/s)
        HCNalt: float, optional
            Altitude HCN is produced in the atmosphere (cm)
        nz : integer, optional
            number of vertical layers in the atmosphere for HCN diffusion
            calculation
        T_trop : float, optional
            Tropopause temperature (K)
        P_trop: float, optional
            Tropopause pressure (K)
        vo : float, optional
            Turnover velocity of the ocean (cm/s)
        zs : float, optional
            Depth of the surface ocean (cm)
        zd : float, optional
            Depth of the deep ocean (cm)
        fH2O_trop: type
            H2O mixing ratio at the tropopause
        """
        self.species = ['H2O','H2','CO','CO2','CH4','N2','NH3', \
                        'HCN', 'C2Hn', 'Haze','HCONH2','HCOOH']
        self.Navo = 6.022e23

        # atmosphere parameters
        self.T_surf = T_surf # Surface T (K)
        # ^ Here we assume the surface temperature of
        # atmosphere is the same as ocean temperature (seems ok)
        self.pH_ocean = pH_ocean # ocean pH (no units)

        # optional parameters
        self.Kzz = Kzz # eddy diffusion (cm2/s)
        self.HCNalt = HCNalt # altitude HCN is produced (cm)
        self.nz = nz # number of vertical layers in atmosphere
        self.T_trop = T_trop # Tropopause temperature (K)
        atmos.t_trop = T_trop
        self.P_trop = P_trop # Tropopause pressure (bar)
        atmos.p_trop = P_trop
        self.vo = vo # Turnover velocity of ocean (cm/s)
        self.zs = zs # depth of surface ocean (cm)
        self.zd = zd # depth of deep ocean (cm)
        atmos.fh2o = fH2O_trop # tropopause H2O mixing ratio

        # Lake parameters
        self.T_lake = T_lake # lake temperature (K)
        self.pH_lake = pH_lake # lake pH (no units)
        self.ktot_lake = HCN_hydrolysis_rate(T_lake, pH_lake) # 1/s
        self.kform_lake = HCONH2_hydrolysis_rate(T_lake,pH_lake) # 1/s

        # HCN henry's law constant
        lnkH = 8205.7/self.T_lake-25.323 # in M/atm
        self.kH_lake = np.exp(lnkH)/1.013 # conver to M/bar

        # equilibrium constant
        self.pKa_HCN_lake = -8.85+3802./self.T_lake+.01786*self.T_lake

    def rhs_verbose(self,t,y):
        N = y[:-2]
        HCONH2, HCOOH = y[-2:]/self.Navo

        # the right-hand-side of the atmosphere
        dNdt, pressure, T_surf, tau_uv, N_H2O, mubar = atmos.rhs_verbose(t,N)

        # compute diffusion of HCN to the surface
        alt, fHCN = diffuseHCN(dNdt[self.species.index('HCN')], Ts = self.T_surf, \
                            Ps = pressure, mubar = mubar, pH = self.pH_ocean, \
                            Kzz = self.Kzz, top_atm = self.HCNalt, nz = self.nz, \
                            T_trop = self.T_trop, P_trop = self.P_trop, vo = self.vo, \
                            zs = self.zs, zd = self.zd)

        # surface partial pressure of HCN
        PHCN = fHCN[0]*pressure

        # dissolved HCN in the lake on a little island
        HCN_aq = self.kH_lake*PHCN

        # dissolved CN in the lake
        CN = HCN_aq*10.**(self.pH_lake-self.pKa_HCN_lake)

        # total HCN in the lake
        sumHCN = HCN_aq + CN

        # hydrolysis of HCN and HCONH2 in the lake
        dHCONH2_dt = (self.ktot_lake*sumHCN - self.kform_lake*HCONH2)*self.Navo
        dHCOOH_dt = (self.kform_lake*HCONH2)*self.Navo

        # right-hand-side of ODEs
        RHS = np.append(dNdt,np.array([dHCONH2_dt,dHCOOH_dt]))

        return RHS, pressure, T_surf, tau_uv, N_H2O, mubar, PHCN, HCN_aq, CN, sumHCN

    def rhs(self,t,y):
        RHS, pressure, T_surf, tau_uv, N_H2O, mubar, PHCN, HCN_aq, CN, sumHCN = self.rhs_verbose(t,y)
        return RHS

    def integrate(self,init_dict,tspan=[0,np.inf],H2end=1e-2,method = "LSODA",rtol=1e-6,**kwargs):
        '''Evolves a Hadean Earth atmosphere using a simple 0-D photochemical
        model coupled to a surface lake from time tspan[0] to tspan[1]
        given the initial conditions Ninit_dict.

        Parameters
        ----------
        tspan : list
            Will integrate from time tspan[0] to tspan[1] in seconds.
        Ninit_dict : dict
            Dictionary containing initial conditions. Must contain
            the following 11 dictionary items: H2, CO, CO2, CH4, N2, NH3, HCN, C2Hn,
            Haze, HCONH2, HCOOH. The species HCONH2, HCOOH are mol/kg, in the lake.
            The rest of the species are in molecules/cm2 in the atmosphere.
        out_dict: bool
            If true, then the output will be a dictionary. If false, then the output
            will be a numpy array with the solution.
        method: str, optional
            Method used by the scipy integrator.
        **kwargs:
            The same optional arguments as scipy.integrate.solve_ivp,
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

        Returns
        -------
        out: dict
            Dictionary defining the state of the atmosphere and lake as a function of time.
        '''
        init = np.zeros(len(self.species))
        init[0] = 0. # H2O doesn't matter
        init[1:-5] = np.array([init_dict[key]*self.Navo for key in self.species[1:-5]])
        # init[-2] = init_dict['HCONH2']*self.Navo
        # init[-1] = init_dict['HCOOH']*self.Navo

        # stop integration when H2 escapes
        def event(t,y):
            return y[1]-H2end*self.Navo
        event.terminal = True

        # set tau uv to 100. then integrate
        atmos.tau_uv_init = 100.
        sol = solve_ivp(self.rhs,tspan,init,method = method, events=event, rtol=rtol,**kwargs)

        # check if solver succeeded
        if not sol.success:
            sys.exit('ODE integration failed')

        # only take positive solution (it should be almost positive)
        yvals = np.clip(sol.y.T,0,np.inf)
        tvals = sol.t

        # re-run the RHS of the ODEs but this time save interesting stuff
        dydt = np.zeros(yvals.shape)
        pressure = np.zeros(len(tvals))
        T_surf = np.zeros(len(tvals))
        tau_uv = np.zeros(len(tvals))
        N_H2O = np.zeros(len(tvals))
        mubar = np.zeros(len(tvals))
        PHCN = np.zeros(len(tvals))
        HCN_aq = np.zeros(len(tvals))
        CN = np.zeros(len(tvals))
        sumHCN = np.zeros(len(tvals))
        atmos.tau_uv_init = 100.
        for i in range(0,len(tvals)):
            dydt[i],pressure[i],T_surf[i],tau_uv[i],N_H2O[i],mubar[i], PHCN[i], \
            HCN_aq[i], CN[i], sumHCN[i] = \
            self.rhs_verbose(tvals[i],yvals[i])

        # build a dictionary of output
        yvals_mol = yvals/self.Navo
        out = {}
        out['Ntot'] = np.sum(yvals_mol[:,:-5],axis=1)-\
                      yvals_mol[:,self.species.index('H2O')] + N_H2O/self.Navo
        out['Psurf'] = pressure
        out['Tsurf'] = T_surf
        out['tau_uv'] = tau_uv
        out['mubar'] = mubar
        out['time'] = tvals
        out['dNHCN_dt'] = dydt[:,self.species.index('HCN')]
        out['dNHaze_dt'] = dydt[:,self.species.index('Haze')]
        out['NH2O_strat'] = yvals_mol[:,self.species.index('H2O')]
        out['H2O'] = N_H2O/self.Navo/out['Ntot']
        for i in range(1,len(self.species)-5):
            out[self.species[i]] = yvals_mol[:,i]/out['Ntot']
        out['sumHCN'] = yvals_mol[:,-5]
        out['sumC2Hn'] = yvals_mol[:,-4]
        out['sumHaze'] = yvals_mol[:,-3]
        out['PHCN'] = PHCN
        out['HCN_aq_lake'] = HCN_aq
        out['CN_lake'] = CN
        out['sumHCN_lake'] = sumHCN
        out['HCONH2_lake'] = yvals[:,-2]/self.Navo
        out['HCOOH_lake'] = yvals[:,-1]/self.Navo

        return out

    def loadgas(self,gas):
        self.stm = SteamAtm(gas)
        self.stm.verbose = False

    def impact_integrate(self,N_H2O_ocean,N_CO2,N_N2,M_i,**kwargs):
        solution = self.stm.impact(N_H2O_ocean,N_CO2,N_N2,M_i)
        Ninit_dict = {}
        for i in range(1,len(self.species)-5):
            Ninit_dict[self.species[i]] = solution[self.species[i]][-1]*solution['Ntot'][-1]
        Ninit_dict['NH3'] = 0.0
        out = self.integrate(Ninit_dict,**kwargs)
        return out
