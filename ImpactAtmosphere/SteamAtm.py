import numpy as np
import cantera as ct

from .SteamAtmBase import SteamAtmBase
from . import constants as const
from .utils import sat_pressure_H2O, OLR

class SteamAtm(SteamAtmBase):

    def __init__(self, gas, T_prime = 2000, impactor_energy_frac = 0.5, \
        Fe_react_frac = 1, impactor_Fe_frac = 0.33, v_i = 17e5, 
        Ni_area = 1.0):

        # Initialize base class
        SteamAtmBase.__init__(self, gas, T_prime, impactor_energy_frac, 
                             Fe_react_frac, impactor_Fe_frac, v_i, Ni_area)

        # Settings specific to SteamAtm
        self.rtol = 1e-10 # relative tolerance for CVODES integrator
        self.atol = 1e-20 # absolute tolerance for CVODES integrator
        self.dTemp = 3 # Temperature descritization (K)
        self.T_stop = 400
        self.verbose = True 

    def impact(self, N_H2O_ocean, N_CO2, N_N2, M_i, N_CO = 0.0, N_H2 = 0.0, N_CH4 = 0.0, include_condensing_phase = True):
        """Simulates chemistry of a cooling steamy atmosphere after large asteroid impact

        Parameters
        ----------
        N_H2O_ocean : float
            Initial column abundance of H2O (mol/cm2)
        N_CO2 : float
            Initial column abundance of CO2 (mol/cm2)
        N_N2 : float
            Initial column abundance of N2 (mol/cm2)
        M_i : type
            Mass of impactor (g)
        N_CO : float
            Initial column abundance of CO (mol/cm2)
        N_H2 : float
            Initial column abundance of H2 (mol/cm2)
        N_CH4 : float
            Initial column abundance of CH4 (mol/cm2)

        Returns
        -------
        solution : dict
            Dictionary containing the composition of the atmosphere as a function
            of time.
        """

        N_init, P_init, X = \
            self.initial_conditions(N_H2O_ocean,N_CO2,N_N2,M_i,N_CO, N_H2, N_CH4)
        sol = self.cooling_steam_atmosphere_1(N_init)
        if include_condensing_phase:
            sol = self.cooling_steam_atmosphere_2(sol)
        return sol.to_dict()

    def cooling_steam_atmosphere_1(self, N_init):
        """Integrates atmosphere from a very hot state (2000 K) down to where 
        H2O begins to condense.
        """
        
        # create solution object and append initial conditions
        sol = SteamAtmSolution(self)
        tn = 0.0
        Tsurf_prev = self.T_prime
        Tsurf = self.T_prime
        N_prev = N_init.copy()
        N = N_init.copy()
        sol.append(tn, Tsurf, N)

        while True:
            if self.verbose:
                print('Atmospheric Temperature = '+"{:<8}".format('%.1f'%Tsurf_prev+' K'),end='\r')
            
            # Atmospheric properties at the beginning of the big time step
            mubar_prev, Psurf_prev, Ntot_prev, Mtot_prev, mix_prev = self.prep_atmosphere(Tsurf_prev, N_prev)
            dtime = self.dtime_dTemp_first(Tsurf_prev, N_prev) # get the size of the time-step
            self.gas.TPX = Tsurf_prev,Psurf_prev/10.0,mix_prev
            
            # advance chemistry forward dtime using cantera
            r = ct.ConstPressureReactor(self.gas,energy='off')
            if self.surface_catalyst:
                self.surf_phase.TPX = Tsurf_prev,Psurf_prev/10.0,{'Ni':1}
                # advance material surface for 10 seconds while keeping
                # gas phase fixed. This makes the material properties
                # realistic and consistent with the gas
                self.surf_phase.advance_coverages(1.0e10,max_steps=1e7,rtol=self.rtol,atol=self.atol)
                # volume of a 1 cm2 column of the atmosphere. I multiply scale height
                # by 1.0 cm2
                Hscale_prev = const.N_avo*const.k_boltz*Tsurf_prev/(mubar_prev*self.grav)
                volume_column = Hscale_prev*1.0
                r.volume = volume_column/1.0e6 # convert to m3
                Ni_area_SI = self.Ni_area/1.0e4 # conver to m2
                # Make the catalytic nickel surface part of the reactor network
                rs = ct.ReactorSurface(kin=self.surf_phase, r=r, A=Ni_area_SI)
            reactorNetwork = ct.ReactorNet([r])
            reactorNetwork.rtol = self.rtol
            reactorNetwork.atol = self.atol
            reactorNetwork.max_steps = 1000000
            reactorNetwork.advance(dtime)
            self.gas.X = self.gas.X.clip(min=0) # force positive solution

            # save the solution
            tn += dtime
            Tsurf -= self.dTemp
            mubar = self.gas.mean_molecular_weight
            Ntot = Mtot_prev/mubar
            N[:] = Ntot*self.gas.X.clip(min=0)
            # break the integration if H2O is supersaturated
            if Tsurf < const.T_crit_H2O:
                mubar, Psurf, Ntot, Mtot, mix = self.prep_atmosphere(Tsurf, N)
                P_H2O_sat = sat_pressure_H2O(Tsurf)
                P_H2O = Psurf*mix[self.ind_H2O]
                if P_H2O > P_H2O_sat:
                    # Adjust N_H2O so that it is at saturation
                    N_H2O = self.saturated_H2O_column(Tsurf, N)
                    N[self.ind_H2O] = N_H2O
                    # we append the adjusted column
                    sol.append(tn, Tsurf, N)
                    break 
            
            sol.append(tn, Tsurf, N) # save

            # update _prev
            Tsurf_prev = Tsurf
            N_prev[:] = N[:]

        return sol

    def cooling_steam_atmosphere_2(self, sol):
        """Evolves a cooling steam atmosphere where H2O is condensing. 
        `sol` gives the initial conditions
        """

        # initial conditions
        tn = sol.time[-1]
        Tsurf_prev = sol.Tsurf[-1]
        Tsurf = sol.Tsurf[-1]
        N = np.empty(self.gas.n_total_species)
        for j,name in enumerate(self.gas.species_names):
            N[j] = sol.mix[name][-1]*sol.Ntot[-1]
        N_prev = N.copy()

        while True:
            if self.verbose:
                print('Atmospheric Temperature = '+"{:<8}".format('%.1f'%Tsurf_prev+' K'),end='\r')

            # Atmospheric properties at the beginning of the big time step
            mubar_prev, Psurf_prev, Ntot_prev, Mtot_prev, mix_prev = self.prep_atmosphere(Tsurf_prev, N_prev)
            dtime = self.dtime_dTemp_second(Tsurf_prev, N_prev) # get the size of the time-step
            self.gas.TPX = Tsurf_prev,Psurf_prev/10.0,mix_prev

            # advance chemistry forward dtime using cantera
            r = ct.ConstPressureReactor(self.gas,energy='off')
            if self.surface_catalyst:
                self.surf_phase.TPX = Tsurf_prev,Psurf_prev/10.0,{'Ni':1}
                # advance material surface for 10 seconds while keeping
                # gas phase fixed. This makes the material properties
                # realistic and consistent with the gas
                self.surf_phase.advance_coverages(1.0e5,max_steps=1e7,rtol=self.rtol,atol=self.atol)
                # volume of a 1 cm2 column of the atmosphere. I multiply scale height
                # by 1.0 cm2
                Hscale_prev = const.N_avo*const.k_boltz*Tsurf_prev/(mubar_prev*self.grav)
                volume_column = Hscale_prev*1.0
                r.volume = volume_column/1.0e6 # convert to m3
                Ni_area_SI = self.Ni_area/1.0e4 # conver to m2
                # Make the catalytic nickel surface part of the reactor network
                rs = ct.ReactorSurface(kin=self.surf_phase, r=r, A=Ni_area_SI)
            reactorNetwork = ct.ReactorNet([r])
            reactorNetwork.rtol = self.rtol
            reactorNetwork.atol = self.atol
            reactorNetwork.max_steps = 1000000
            reactorNetwork.advance(dtime)
            self.gas.X = self.gas.X.clip(min=0) # force positive solution

            # save the solution
            tn += dtime
            Tsurf -= self.dTemp
            mubar = self.gas.mean_molecular_weight
            Ntot = Mtot_prev/mubar
            N[:] = Ntot*self.gas.X.clip(min=0)
            # Adjust N_H2O so that it is at saturation
            N_H2O = self.saturated_H2O_column(Tsurf, N)
            N[self.ind_H2O] = N_H2O
            sol.append(tn, Tsurf, N) # save

            # update _prev
            Tsurf_prev = Tsurf
            N_prev[:] = N[:]

            if Tsurf < self.T_stop:
                break

        return sol

    def saturated_H2O_column(self, T, N):
        """Computes the H2O column abundances, assuming H2O is at
        saturation."""
        P_H2O = sat_pressure_H2O(T)

        N_tmp = N.copy()
        N_tmp[self.ind_H2O] = 0.0
        N_dry = np.sum(N_tmp)
        X_dry = N_tmp/N_dry
        self.gas.X = X_dry
        mu_dry = self.gas.mean_molecular_weight

        PP = np.empty(3)
        PP[0] = const.mu_H2O*self.grav
        PP[1] = N_dry*mu_dry*self.grav - P_H2O
        PP[2] = -N_dry*P_H2O

        N_H2O = np.max(np.real(np.roots(PP)))
        return N_H2O

    def prep_atmosphere(self, Tsurf, N):
        """Computes atmospheric properties from the column abundance
        """
        Ntot = np.sum(N)
        mix = N/Ntot
        self.gas.X = mix
        mubar = self.gas.mean_molecular_weight
        Mtot = Ntot*mubar
        Psurf = Mtot*self.grav
        return mubar, Psurf, Ntot, Mtot, mix

    def dtime_dTemp_first(self, T, N):  
        mubar, Psurf, Ntot, Mtot, mix = self.prep_atmosphere(T, N)
        Fir = OLR(T)
        dT_dt = -self.grav/(const.cp_H2O*Psurf)*Fir
        dtime = -(1/dT_dt)*self.dTemp
        return dtime

    def dtime_dTemp_second(self, T, N):  
        mubar, Psurf, Ntot, Mtot, mix = self.prep_atmosphere(T, N)
        Fir = OLR(T)
        P_H2O = sat_pressure_H2O(T)
        dT_dt = (-self.grav/(const.cp_H2O*Psurf))*Fir \
                *(1 + const.L_H2O**2*const.mu_H2O*P_H2O/(const.cp_H2O*Psurf*const.R*T**2))**-1
        dtime = -(1/dT_dt)*self.dTemp
        return dtime

class SteamAtmSolution():
    """helper class to save results
    """

    def __init__(self, stm):
        self.stm = stm
        self.time = np.empty((0,),np.float64)
        self.Tsurf = np.empty((0,),np.float64)
        self.mubar = np.empty((0,),np.float64)
        self.Psurf = np.empty((0,),np.float64)
        self.Ntot = np.empty((0,),np.float64)
        self.mix = {}
        for i,name in enumerate(self.stm.gas.species_names):
            self.mix[name] = np.empty((0,),np.float64)

    def append(self, t, Tsurf, N):
        mubar, Psurf, Ntot, Mtot, mix = self.stm.prep_atmosphere(Tsurf, N)
        # append
        self.time = np.append(self.time, t)
        self.Tsurf = np.append(self.Tsurf, Tsurf)
        self.mubar = np.append(self.mubar, mubar)
        self.Psurf = np.append(self.Psurf, Psurf)
        self.Ntot = np.append(self.Ntot, Ntot)
        for i,name in enumerate(self.stm.gas.species_names):
            self.mix[name] = np.append(self.mix[name],mix[i])

    def to_dict(self):
        sol = {}
        sol['time'] = self.time
        sol['Tsurf'] = self.Tsurf
        sol['mubar'] = self.mubar
        sol['Psurf'] = self.Psurf
        sol['Ntot'] = self.Ntot
        for key in self.mix:
            sol[key] = self.mix[key]
        return sol