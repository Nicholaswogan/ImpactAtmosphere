import numpy as np
import cantera as ct
from cantera import CanteraError

from .SteamAtmBase import SteamAtmBase
from ._cvode import CVode
from . import constants as const
from . import utils

class SteamAtmContinuous(SteamAtmBase):

    def __init__(self, gas, T_prime = 2000, impactor_energy_frac = 0.5, \
        Fe_react_frac = 1, impactor_Fe_frac = 0.33, v_i = 20.7e5, 
        Ni_area = 1.0):

        # Initialize base class
        SteamAtmBase.__init__(self, gas, T_prime, impactor_energy_frac, 
                             Fe_react_frac, impactor_Fe_frac, v_i, Ni_area)

        # Indexes for melt-atmosphere equilibration
        self.nmelt_1 = 1 + self.ngas + self.nsurf
        self.ind_O2 = self.gas.species_names.index('O2')

        # Settings specific to SteamAtmContinuous
        self.rtol = 1.0e-4 # relative tolerance
        self.atol = 1.0e-20 # absolute tolerance
        self.dTime = 10.0 # frequency in which solution is saved (years)
        self.max_integrator_errors = 10
        self.T_stop = 400.0 # temperature to stop the integration
        self.initial_integration = False

    def impact(self, N_H2O_ocean, N_CO2, N_N2, M_i, N_CO=0.0, N_H2=0.0, N_CH4=0.0, include_condensing_phase=True):
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
        self.initial_conditions(N_H2O_ocean, N_CO2, N_N2, M_i, N_CO, N_H2, N_CH4)
        sol = self.cooling_steam_atmosphere_1(N_init, P_init, X)
        if include_condensing_phase:
            sol = self.cooling_steam_atmosphere_2(sol)
        return sol.to_dict()

    def impact_from_dict(self, init_cond, include_condensing_phase=True):

        N_init = np.empty(len(self.gas.species_names))
        X = np.empty(len(self.gas.species_names))
        for i,name in enumerate(self.gas.species_names):
            N_init[i] = init_cond['Ntot']*init_cond[name]
            X[i] = init_cond[name]
        P_init = init_cond['Psurf']

        sol = self.cooling_steam_atmosphere_1(N_init, P_init, X)
        if include_condensing_phase:
            sol = self.cooling_steam_atmosphere_2(sol)
        return sol.to_dict()

    def equilibrate_initial_conditions(self, N_H2O_ocean, N_CO2, N_N2, M_i, N_CO = 0.0, N_H2 = 0.0, N_CH4 = 0.0,
                                       M_melt = None, DQFM_melt = None, equilibrium_time = 1000.0):
        """Integrates the Atmosphere for a period of time to equilibrate it. The main purpose of this function is to
        allow for equilibration between the atmophere and a melt pond.

        Parameters
        ----------
        N_H2O_ocean : float
            Column abundance of the ocean (mol/cm^2)
        N_CO2 : float
            Column abundance of CO2 in the atmosphere (mol/cm^2)
        N_N2 : float
            Column abundance of N2 in the atmosphere (mol/cm^2)
        M_i : _type_
            Mass of the impactor (grams)
        N_CO : float, optional
            Column abundance of CO in the atmosphere (mol/cm^2), by default 0.0
        N_H2 : float, optional
            Column abundance of H2 in the atmosphere (mol/cm^2), by default 0.0
        N_CH4 : float, optional
            Column abundance of CH4 in the atmosphere (mol/cm^2), by default 0.0
        M_melt : float, optional
            Mass of the melt pond in grams, by default None. If None, no melt pond is assumed.
        DQFM_melt : float, optional
            Redox state of the melt pond relative to the QFM buffer, by default None
        equilibrium_time : float, optional
            Time in which to assume equilibration has be achieved (years), by default 1000.0.

        Returns
        -------
        dict
            A dictionary containing the equilibrated state of the atmosphere.
        """        

        if M_melt != None:
            melt_reactions = True
            if DQFM_melt == None:
                raise Exception("DQFM_melt must be specified if M_melt is specified")
        else:
            melt_reactions = False
            
        N_init, P_init, X = \
        self.initial_conditions(N_H2O_ocean, N_CO2, N_N2, M_i, N_CO, N_H2, N_CH4)
        if melt_reactions:
            mols_melt = utils.alter_basalt_melt_to_fO2(self.T_prime, P_init/1e6, DQFM_melt)
        else:
            mols_melt = None

        if self.surface_catalyst:
            X_surf = self.equilibrium_coverages(self.T_prime, P_init, X)
            y0_surf = X_surf*self.Ni_area*self.site_density*1.0e12 # pico mol catalyst/cm2 Earth
        else:
            y0_surf = np.array([])

        if melt_reactions:
            y0_melt = np.array([mols_melt[0],mols_melt[1], 0.0]) # start with no H2O
        else:
            y0_melt = np.array([])
            
        y0 = np.concatenate(([self.T_prime], N_init, y0_surf, y0_melt))

        args = (self, melt_reactions, mols_melt, M_melt)
        t0 = 0.0
        c = CVode(rhs_equilibrate, t0, y0, rtol=self.rtol, atol=self.atol, args=args, mxsteps=10000)
        t = np.array(0.0)
        y = np.empty(y0.shape[0])
        tn = equilibrium_time*const.yr

        tries = 0
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
            else:
                break

        if ret < 0:
            raise Exception('re_equilibrate_initial_conditions failed: '+str(ret))

        y = np.clip(y, 1.0e-30, np.inf)
        Tsurf, N, mubar, Psurf, Ntot, mix = self.prep_atm_first(y)

        out = {}
        out['Tsurf'] = Tsurf
        out['mubar'] = mubar
        out['Psurf'] = Psurf
        out['Ntot'] = Ntot
        for i,name in enumerate(self.gas.species_names):
            out[name] = mix[i]
        if melt_reactions:
            out['Fe2O3'] = y[self.nmelt_1]
            out['FeO'] = y[self.nmelt_1+1]
            out['H2O_melt'] = y[self.nmelt_1+2]
            fO2 = out['O2']*out['Psurf']/1e6
            out['DFMQ'] = np.log10(fO2) - utils.log10fO2_FMQ(Tsurf) 
        
        return out

    def cooling_steam_atmosphere_1(self, N_init, P_init, X):
        """This routine integrates the atmosphere from T_prime, which is assumed to
        be very hot (~2000 K), until water begins to condense.
        """

        # initial conditions
        t0 = 0.0
        if self.surface_catalyst:
            X_surf = self.equilibrium_coverages(self.T_prime, P_init, X)
            y0_surf = X_surf*self.Ni_area*self.site_density*1.0e12 # pico mol catalyst/cm2 Earth
        else:
            y0_surf = np.array([])
            
        y0 = np.concatenate(([self.T_prime], N_init, y0_surf))

        # integrate initial conditions at constant T for some time.
        # This is needed for when there is a redox buffer, which
        # is not accounted for during chem equil calculation.
        # if self.redox_buffer:
        #     y0 = self.re_equilibrate_initial_conditions(y0)
        
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

    def re_equilibrate_initial_conditions(self, y0):
        """The purpose of this routine is to integrate the inital state for some time to
        to make sure it is at equilibrium. This is only necessary for cases that use
        a redox buffer (e.g. FMQ)
        """
        self.initial_integration = True
        args = (self,)
        t0 = 0.0
        c = CVode(rhs_first, t0, y0, rtol=self.rtol, atol=self.atol, args=args, mxsteps=10000)
        t = np.array(0.0)
        y = np.empty(y0.shape[0])
        tn = 1000.0*const.yr # 1000 years

        tries = 0
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
            else:
                break

        if ret < 0:
            raise Exception('re_equilibrate_initial_conditions failed: '+str(ret))
        y = np.clip(y, 1.0e-30, np.inf)
        self.initial_integration = False
        return y
    
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
        N = y[self.ngas_1:self.nsurf_1]
        Ntot = np.sum(N)
        mix = N/Ntot
        self.gas.X = mix
        mubar = self.gas.mean_molecular_weight
        Psurf = Ntot*mubar*self.grav
        return T, N, mubar, Psurf, Ntot, mix

    def prep_atm_first_catalyst(self, y):
        N_cat = y[self.nsurf_1:self.nsurf_1+self.nsurf]
        X_cat = N_cat/np.sum(N_cat)
        return X_cat

    def prep_atm_second(self, y):
        T = y[0]
        P_H2O = utils.sat_pressure_H2O(T)
        
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

def rhs_equilibrate(t, y, dy, self, melt_reactions, mols_melt, M_melt):
    """Right-hand-side function for equilibrating the atmosphere with a pond of magma. Here I implement
    to fake reactions: 
    0.5O2 + 2 FeO <=> Fe2O3
    H2O(g) <=> H2O(aq)
    The first reaction pushes and pulls O2 out of the magma so that there is redox-equilibrium between
    the melt and atmosphere. The equilibrium is defined by a relation in Kress and Carmichael (1991).
    The second reaction pushes and pulls H2O out of the magma to satisfy the H2O solubility in magma.
    In both cases, the reaction rates are made up, and are designed to be fast to enforce an equilibrium
    state on a reasonable timescale (seconds).
    """

    T, N, mubar, Psurf, Ntot, mix = self.prep_atm_first(y)
    
    # This computes contributions from chemistry (gas-phase and surfaces)
    ret = rhs_first(t, y, dy, self)

    # Add rates of change from reactions with the melt
    if melt_reactions:
        n_Fe2O3, n_FeO, n_H2O_melt = y[self.nmelt_1:]
        fO2 = (Psurf/1e6)*mix[self.ind_O2] # Oxygen pressure in bars

        # reaction driving redox between atmosphere and melt to equilibrium
        dn_O2_dt, dn_FeO_dt, dn_Fe2O3_dt = utils.rates_melt_reactions(T, Psurf/1e6, n_Fe2O3, n_FeO, fO2, mols_melt)
        dy[1+self.ind_O2] = dy[1+self.ind_O2] + dn_O2_dt*M_melt/self.area # mol/cm^2 Earth
        dy[self.nmelt_1] = dn_Fe2O3_dt
        dy[self.nmelt_1+1] = dn_FeO_dt

        # reaction driving H2O soubility in melt to equilibrium
        P_H2O = mix[self.ind_H2O]*Psurf/1e6 # bars
        dn_H2O_gas_dt, dn_H2O_melt_dt = utils.rates_melt_H2O_reaction(P_H2O, n_H2O_melt)
        dy[1+self.ind_H2O] = dy[1+self.ind_H2O] + dn_H2O_gas_dt*M_melt/self.area # mol/cm^2 Earth
        dy[self.nmelt_1+2] = dn_H2O_melt_dt

    # temperature does not change
    dy[0] = 0.0

    return ret

def rhs_first(t, y, dy, self):
    """right-hand-side function for when the atmosphere is
    is very hot, and has no condensing H2O.
    """   
    dy[:] = 0.0  

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
        dy[self.nsurf_1:self.nmelt_1] = surf_rates*self.Ni_area*1.0e12 # pico mol/(cm2 Earth * s)
        dy[self.ngas_1:self.nsurf_1] = gas_rates*self.Ni_area # mol/(cm2 Earth * s)

    # gas-phase rates
    gas_rates = self.gas.net_production_rates*1.0e3/1.0e6 # mol/(cm3 * s)
    dy[self.ngas_1:self.nsurf_1] = dy[self.ngas_1:self.nsurf_1] + gas_rates[:]*Ha # moles/cm2/s

    # climate
    Fir = utils.net_outgoing_flux(T)
    dT_dt = -self.grav/(const.cp_H2O*Psurf)*Fir
    if self.initial_integration:
        # we don't evolve temperature
        dy[0] = 0.0
    else:
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
        gout[0] = utils.sat_pressure_H2O(T) - PH2O
    else:
        gout[0] = 1.0

    # for some reason integration fails at 1500 K
    # So I stop and restart beyond 1500 K
    gout[1] = T - (1500.0 + const.TOL)
    return 0

def rhs_second(t, y, dy, self):
    """right-hand-side function for when the atmosphere has condensing H2O.
    """  
    dy[:] = 0.0 

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
    
    # gas rates
    gas_rates = self.gas.net_production_rates*1.0e3/1.0e6 # mol/(cm3 * s)
    dN_dt = gas_rates*Ha # moles/cm2/s
    dy[1:self.ind_H2O+1] = dy[1:self.ind_H2O+1] + dN_dt[:self.ind_H2O]
    dy[self.ind_H2O+1:self.ngas] = dy[self.ind_H2O+1:self.ngas] + dN_dt[self.ind_H2O+1:]
    
    # climate
    Fir = utils.net_outgoing_flux(T)
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