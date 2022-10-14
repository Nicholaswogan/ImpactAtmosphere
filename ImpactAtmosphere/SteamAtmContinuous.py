import numpy as np
import cantera as ct
from cantera import CanteraError

from .SteamAtmBase import SteamAtmBase
from ._cvode import CVode
from . import constants as const
from .utils import sat_pressure_H2O, OLR

class SteamAtmContinuous(SteamAtmBase):

    def __init__(self, gas, T_prime = 2000, impactor_energy_frac = 0.5, \
        Fe_react_frac = 1, impactor_Fe_frac = 0.33, v_i = 17e5, 
        Ni_area = 1.0):

        # Initialize base class
        SteamAtmBase.__init__(self, gas, T_prime, impactor_energy_frac, 
                             Fe_react_frac, impactor_Fe_frac, v_i, Ni_area)

        # Settings specific to SteamAtmContinuous
        self.rtol = 1.0e-4 # relative tolerance
        self.atol = 1.0e-20 # absolute tolerance
        self.dTime = 10.0 # frequency in which solution is saved (years)
        self.max_integrator_errors = 10
        self.T_stop = 400.0

    def impact(self,N_H2O_ocean,N_CO2,N_N2,M_i,N_CO = 0.0, N_H2 = 0.0, N_CH4 = 0.0, include_condensing_phase = True):
        """Simulates chemistry of a cooling steamy atmosphere after large asteroid
        impact. The simulation stops when the steam begins to condense.

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
        sol = self.cooling_steam_atmosphere_1(N_H2O_ocean, N_init, P_init, X)
        if include_condensing_phase:
            sol = self.cooling_steam_atmosphere_2(sol)
        return sol.to_dict()

    def cooling_steam_atmosphere_1(self, N_H2O_ocean, N_init, P_init, X):
        """This routine integrates the atmosphere from T_prime, which is assumed to
        be very hot (~2000 K), until water begins to condense.
        """

        # initial conditions
        t0 = 0.0
        if self.surface_catalyst:
            X_surf = self.equilibrium_coverages(self.T_prime, P_init, X)
            y0_surf = X_surf*self.Ni_area*self.site_density*1.0e12 # pico mol catalyst/cm2 Earth
            y0 = np.concatenate(([self.T_prime], N_init, y0_surf))
        else:
            y0 = np.append(self.T_prime, N_init)
        
        # initialize integrator
        args = (self,)
        ng = 2
        c = CVode(rhs_first, t0, y0, rtol=self.rtol, atol=self.atol, args=args, ng=ng, g=stop_first)
        
        # set up ability to save solution
        tn = 0.0
        t = np.array(0.0)
        y = np.empty(y0.shape[0])
        sol = SteamAtmContinuousSolution(self) # solution object
        sol.append(t0, y0, 1) # load initial conditions
        roots = np.zeros(ng,np.int32) # for root finding

        while True:
            tries = 0
            tn += self.dTime*const.yr # advance time by dTime
            while True:
                t[...] = tn
                ret = c.integrate(t, y)
                if ret < 0:
                    # error was encountered.
                    # we re-initialize integrator and try again.
                    if tries > self.max_integrator_errors:
                        # give up
                        break
                    # clip y for good luck
                    y = np.clip(y, 1.0e-30, np.inf)
                    c.reinit(t, y) # re-initialize integrator
                    tries += 1
                elif ret == 0 or ret == 1:
                    # successful return
                    # append solution
                    sol.append(t.item(), y, 1)
                    break
                elif ret == 2:
                    c.rootinfo(roots)
                    if roots[0] != 0:
                        # append the transition
                        sol.append(t.item(), y, 1)
                        break
                    elif roots[1] != 0:
                        # Here we skip over exactly T = 1500
                        # for some reason it causes problems
                        y[0] = 1500.0 - const.TOL
                        c.reinit(t, y)

            if ret < 0 or roots[0] != 0:
                # stop if there is an error, or if 
                # water has begun condensing
                break
        
        if ret < 0:
            raise Exception('First impact integration failed.')

        return sol

    def cooling_steam_atmosphere_2(self, sol):
        """integrates a steam atmosphere that has actively condensing H2O.
        Initial conditions are given in `sol`, and we append the solution to
        `sol`.
        """

        # Initial conditions
        t0 = sol.time[-1]
        N = np.empty(self.ngas - 1)
        i = 0
        for j,name in enumerate(self.gas.species_names):
            if name != 'H2O':
                N[i] = sol.mix[name][-1]*sol.Ntot[-1]
                i += 1
        if self.surface_catalyst:
            y0 = np.concatenate(([sol.Tsurf[-1]], N, sol.y[1+self.ngas:]))
        else:
            y0 = np.append(sol.Tsurf[-1], N)

        # initialize integrator
        args = (self,)
        ng = 1
        c = CVode(rhs_second, t0, y0, rtol=self.rtol, atol=self.atol, args=args, ng=ng, g=stop_second)

        # setup ability to integrate
        tn = t0
        t = np.array(0.0)
        y = np.empty(y0.shape[0])
        roots = np.zeros(ng,np.int32) # for root finding

        while True:
            tries = 0
            tn += self.dTime*const.yr
            while True:
                t[...] = tn
                ret = c.integrate(t, y)
                if ret < 0:
                    # error was encountered.
                    # we re-initialize integrator and try again.
                    if tries > self.max_integrator_errors:
                        # give up
                        break
                    # clip y for good luck
                    y = np.clip(y, 1.0e-30, np.inf)
                    c.reinit(t, y) # re-initialize integrator
                    tries += 1
                elif ret == 0 or ret == 1:
                    # successful return
                    # append solution
                    sol.append(t.item(), y, 2)
                    break
                elif ret == 2:
                    c.rootinfo(roots)
                    if roots[0] != 0:
                        # append the end state
                        sol.append(t.item(), y, 2)
                        break
            if ret < 0 or roots[0] != 0:
                # stop if there is an error, or if 
                # the integration is finished
                break

        if ret < 0:
            raise Exception('Second impact integration failed.')

        return sol
    
    def equilibrium_coverages(self, T, P, X):
        """Given gas-phase conditions, this routine computes
        Equilibrium coverages for the surface.
        """
        self.gas.TPX = T, P/10.0, X
        self.surf_phase.TPX = T, P/10.0, {'Ni':1.0}
        self.surf_phase.advance_coverages(1.0e10, max_steps=1e7)
        return self.surf_phase.X

    def prep_atm_first(self, y):
        T = y[0]
        N = y[1:self.ngas+1]
        Ntot = np.sum(N)
        mix = N/Ntot
        self.gas.X = mix
        mubar = self.gas.mean_molecular_weight
        Psurf = Ntot*mubar*self.grav
        return T, N, mubar, Psurf, Ntot, mix

    def prep_atm_first_catalyst(self, y):
        N_cat = y[self.ngas+1:]
        X_cat = N_cat/np.sum(N_cat)
        return X_cat

    def prep_atm_second(self, y):
        T = y[0]
        P_H2O = sat_pressure_H2O(T)
        
        N = np.empty(self.ngas)
        N[:self.ind_H2O] = y[1:self.ind_H2O+1]
        N[self.ind_H2O] = 0.0
        N[self.ind_H2O+1:] = y[self.ind_H2O+1:self.ngas]

        # below, we compute the N_H2O
        N_dry = np.sum(N)
        X_dry = N/N_dry
        self.gas.X = X_dry
        mu_dry = self.gas.mean_molecular_weight

        PP = np.empty(3)
        PP[0] = const.mu_H2O*self.grav
        PP[1] = N_dry*mu_dry*self.grav - P_H2O
        PP[2] = -N_dry*P_H2O

        N_H2O = np.max(np.real(np.roots(PP)))
        N[self.ind_H2O] = N_H2O
        Ntot = np.sum(N)
        mix = N/Ntot

        self.gas.X = mix
        mubar = self.gas.mean_molecular_weight
        Psurf = Ntot*mubar*self.grav

        return T, N, P_H2O, mubar, Psurf, Ntot, mix

    def prep_atm_second_catalyst(self, y):
        N_cat = y[self.ngas:]
        X_cat = N_cat/np.sum(N_cat)
        return X_cat

##############################
### RHS and Root functions ###
##############################

# right-hand-side and root functions for
# both integrations involved in the post-impact atmosphere

def rhs_first(t, y, dy, self):
    """right-hand-side function for when the atmosphere is
    is very hot, and has no condensing H2O.
    """    

    try:
        T, N, mubar, Psurf, Ntot, mix = self.prep_atm_first(y)
        self.gas.TPX = T, Psurf/10.0, mix
        if self.surface_catalyst:
            X_cat = self.prep_atm_first_catalyst(y)
            self.surf_phase.TPX = T,Psurf/10.0,X_cat
    except CanteraError:
        return 1
    
    # scale height
    Ha = const.N_avo*const.k_boltz*T/(mubar*self.grav)
    
    if self.surface_catalyst:
        all_rates = self.surf_phase.net_production_rates
        surf_rates = all_rates[:self.nsurf]*1.0e3/1.0e4 # mol/(cm2 metal * s)
        gas_rates = all_rates[self.nsurf:]*1.0e3/1.0e4 # mol/(cm2 metal * s)

        # surface rates
        dy[self.ngas+1:] = surf_rates*self.Ni_area*1.0e12 # pico mol/(cm2 Earth * s)
        dy[1:self.ngas+1] = gas_rates*self.Ni_area # mol/(cm2 Earth * s)
    else:
        gas_rates = self.gas.net_production_rates*1.0e3/1.0e6 # mol/(cm3 * s)
        dy[1:self.ngas+1] = gas_rates*Ha # moles/cm2/s

    # climate
    Fir = OLR(T)
    dT_dt = -self.grav/(const.cp_H2O*Psurf)*Fir
    dy[0] = dT_dt
    
    return 0

def stop_first(t, y, gout, self):
    """root-finding function for when the atmosphere is
    is very hot, and has no condensing H2O. We seek roots
    water begins condensing, and when T = 1500 K, which
    the integrator does not like, so we skip over that exact
    temperature.
    """   

    try:
        T, N, mubar, Psurf, Ntot, mix = self.prep_atm_first(y)
    except CanteraError:
        return 1

    PH2O = Psurf*mix[self.ind_H2O]

    if T < const.T_crit_H2O:
        gout[0] = sat_pressure_H2O(T) - PH2O
    else:
        gout[0] = 1.0

    # for some reason integration fails at 1500 K
    # So I stop and restart beyond 1500 K
    gout[1] = T - (1500.0 + const.TOL)
    return 0

def rhs_second(t, y, dy, self):
    """right-hand-side function for when the atmosphere has condensing H2O.
    """   

    try:
        T, N, P_H2O, mubar, Psurf, Ntot, mix = self.prep_atm_second(y)
        self.gas.TPX = T, Psurf/10.0, mix
        if self.surface_catalyst:
            X_cat = self.prep_atm_second_catalyst(y)
            self.surf_phase.TPX = T,Psurf/10.0,X_cat
    except CanteraError:
        return 1

    # scale height
    Ha = const.N_avo*const.k_boltz*T/(mubar*self.grav)

    if self.surface_catalyst:
        all_rates = self.surf_phase.net_production_rates
        surf_rates = all_rates[:self.nsurf]*1.0e3/1.0e4 # mol/(cm2 metal * s)
        gas_rates = all_rates[self.nsurf:]*1.0e3/1.0e4 # mol/(cm2 metal * s)

        # surface rates
        dy[self.ngas:] = surf_rates*self.Ni_area*1.0e12 # pico mol/(cm2 Earth * s)
        # gas rates
        dN_dt = gas_rates*self.Ni_area # mol/(cm2 Earth * s)
        dy[1:self.ind_H2O+1] = dN_dt[:self.ind_H2O]
        dy[self.ind_H2O+1:self.ngas] = dN_dt[self.ind_H2O+1:]
    else:
        gas_rates = self.gas.net_production_rates*1.0e3/1.0e6 # mol/(cm3 * s)
        dN_dt = gas_rates*Ha # moles/cm2/s
        dy[1:self.ind_H2O+1] = dN_dt[:self.ind_H2O]
        dy[self.ind_H2O+1:self.ngas] = dN_dt[self.ind_H2O+1:]
    
    # climate
    Fir = OLR(T)
    dT_dt = (-self.grav/(const.cp_H2O*Psurf))*Fir \
            *(1 + const.L_H2O**2*const.mu_H2O*P_H2O/(const.cp_H2O*Psurf*const.R*T**2))**-1
    dy[0] = dT_dt
    
    return 0

def stop_second(t, y, gout, self):
    """root-finding function for when the atmosphere has condensing H2O.
    We stop when T == T_stop
    """
    gout[0] = y[0] - self.T_stop
    return 0

class SteamAtmContinuousSolution():
    """A helper class to store the solution as we go along
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

    def append(self, t, y, integ_type):

        if integ_type == 1:
            Tsurf, N, mubar, Psurf, Ntot, mix = self.stm.prep_atm_first(y)
        elif integ_type == 2:
            ind_H2O = self.stm.gas.species_names.index('H2O')
            Tsurf, N, P_H2O, mubar, Psurf, Ntot, mix = self.stm.prep_atm_second(y)
        else:
            raise ValueError('Invalid integ_type')

        # append
        self.time = np.append(self.time, t)
        self.Tsurf = np.append(self.Tsurf, Tsurf)
        self.mubar = np.append(self.mubar, mubar)
        self.Psurf = np.append(self.Psurf, Psurf)
        self.Ntot = np.append(self.Ntot, Ntot)
        for i,name in enumerate(self.stm.gas.species_names):
            self.mix[name] = np.append(self.mix[name],mix[i])

        # we save the current y
        self.y = y

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