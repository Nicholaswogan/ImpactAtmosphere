import numpy as np
import cantera as ct
from cantera import CanteraError
import os

from ._cvode import CVode
from . import constants as const
from .utils import sat_pressure_H2O, OLR

class SteamAtmContinuous():
    def __init__(self, gas, T_prime= 2000, impactor_energy_frac = 0.5, \
        Fe_react_frac = 1, impactor_Fe_frac = 0.33, v_i = 17e5, 
        Ni_area = 1.0):

        self.surface_catalyst = False
        if type(gas) == str:
            zahnle_path = os.path.dirname(os.path.realpath(__file__))+'/data/'
            ct.add_directory(zahnle_path)
            if gas == "Methanation_Ni.yaml":
                # we have a surface catalyst
                self.surface_catalyst = True
                self.surf_phase = ct.Interface('Methanation_Ni.yaml','Ni-surface')
                self.gas = self.surf_phase.adjacent["gas"]
            else:
                self.gas = ct.Solution(gas)
        else:
            self.gas = gas # cantera gas object
        self.gas.basis = 'mass'

        self.ngas = self.gas.n_total_species
        if self.surface_catalyst:
            self.nsurf = self.surf_phase.n_total_species - self.ngas
        else:
            self.nsurf = 0

        # index of H2O
        self.ind_H2O = self.gas.species_names.index('H2O')

        # Parameters
        self.T_prime = T_prime # Initial atmosphere temperature (K)
        self.impactor_energy_frac = impactor_energy_frac # fraction of kinetic energy used to heat the ocean and atmosphere
        self.Fe_react_frac = Fe_react_frac # fraction of the iron that reacts with atmosphere.
        self.impactor_Fe_frac = impactor_Fe_frac # Fe mass fraction of the impactor.
        self.v_i = v_i # velocity of impactor (cm/s)
        if self.surface_catalyst:
            self.Ni_area = Ni_area # cm^2 Ni surface / cm^2 column of the atmosphere
            self.site_density = self.surf_phase.site_density*1.0e3/1.0e4 # mol catalyst/cm2 metal

        # Settings
        self.rtol = 1.0e-4 # relative tolerance
        self.atol = 1.0e-20 # absolute tolerance
        self.dTime = 10.0 # years
        self.max_integrator_errors = 10
        self.T_stop = 400.0

        # Planet
        self.grav = 981.0 # gravity of Earth (cm/s2)
        self.area = 5.1e18 # Area of Earth (cm2)


    def impact(self,N_H2O_ocean,N_CO2,N_N2,M_i,N_CO = 0.0, N_H2 = 0.0, N_CH4 = 0.0):
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

        # P_init = [dynes], N_init = [mol/cm2]
        N_init, P_init, X = \
        self.initial_conditions(N_H2O_ocean,N_CO2,N_N2,M_i,N_CO, N_H2, N_CH4)
        sol = self.cooling_steam_atmosphere_1(N_H2O_ocean, N_init, P_init, X)
        # sol = self.cooling_steam_atmosphere_2(sol)
        return sol.to_dict()

    def cooling_steam_atmosphere_1(self, N_H2O_ocean, N_init, P_init, X):
        """This routine integrates the atmosphere from T_prime, which is assumed to
        be very hot (~2000 K), until water begins to condense.
        """

        # initialize integrator
        t0 = 0.0
        if self.surface_catalyst:
            X_surf = self.equilibrium_coverages(self.T_prime,P_init,X)
            
            # pico mol catalyst/cm2 Earth
            # like the total column
            y0_surf = X_surf*self.Ni_area*self.site_density*1.0e12
            y0 = np.concatenate(([self.T_prime], N_init, y0_surf))
        else:
            y0 = np.append(self.T_prime, N_init)
        args = (self,)
        ng = 2
        c = CVode(rhs_first, t0, y0, rtol=self.rtol, atol=self.atol, args=args, ng=ng, g=stop_first)
        
        # set up ability to save solution
        tn = 0.0
        t = np.array(0.0)
        y = np.empty(y0.shape[0])
        sol = ImpactSolution(self) # solution object
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
        Initial conditions are given in `sol`.
        """

        # initialize integrator
        t0 = sol.time[-1]
        N = np.empty(self.gas.n_total_species - 1)
        i = 0
        for j,name in enumerate(self.gas.species_names):
            if name != 'H2O':
                N[i] = sol.mix[name][-1]*sol.Ntot[-1]
                i += 1
        y0 = np.append(sol.Tsurf[-1], N) # we ignore H2O
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
    
    def equilibrium_coverages(self, T, P, X, time = 100):

        # gas content
        self.gas.TPX = T, P/10.0, X

        ind = self.surf_phase.species_names.index('Ni')
        y0 = np.zeros_like(self.surf_phase.X)
        y0[ind] = 1.0
        self.surf_phase.TPX = T, P/10.0, y0

        self.surf_phase.advance_coverages(time)

        return self.surf_phase.X

    def initial_conditions(self,N_H2O_ocean,N_CO2,N_N2,M_i, N_CO_init, N_H2_init, N_CH4_init):
        """Determines the initial conditions for method `impact`.
        """

        # Determine how much steam is made
        N_H2O_steam = self.steam_from_impact(N_H2O_ocean,N_CO2,N_N2,M_i)
        # React the atmosphere with impactor Fe
        N_H2,N_H2O,N_CO,N_CO2 = self.react_iron(N_H2O_steam,N_CO2,M_i)
        # set as initial conditions
        N_init = np.zeros(len(self.gas.species_names))
        N_init[self.gas.species_names.index('H2')] = N_H2 + N_H2_init
        N_init[self.gas.species_names.index('H2O')] = N_H2O
        N_init[self.gas.species_names.index('CO')] = N_CO + N_CO_init
        N_init[self.gas.species_names.index('CO2')] = N_CO2
        N_init[self.gas.species_names.index('N2')] = N_N2
        N_init[self.gas.species_names.index('CH4')] = N_CH4_init

        P_init = np.sum(self.gas.molecular_weights*N_init*self.grav) # (dynes)
        self.gas.TPX = self.T_prime,P_init*0.1,(N_init/np.sum(N_init))
        self.gas.equilibrate('TP') # equilibrate
        N_init = self.gas.X*np.sum(N_init)
        return N_init, P_init, self.gas.X

    def partial_cp_cgs(self,spec,T):
        '''Finds specific heat capacity of spec in
        ergs/g/K (cgs units).
        '''
        self.gas.TP = T,1e5 # pressure doesn't matter
        index = self.gas.species_names.index(spec)
        cp = self.gas.partial_molar_cp[index]/self.gas.molecular_weights[index]*1e4
        return cp

    def steam_from_impact(self, N_H2O_ocean, N_CO2, N_N2, m_i):
        """Calculates the amount of steam heated to self.T_prime by an impactor
        of mass m_i (grams)

        Parameters
        ----------
        N_H2O_ocean : float
            Size of the ocean (mol/cm2)
        N_CO2 : float
            CO2 in the atmosphere (mol/cm2)
        N_N2 : float
            N2 in the atmosphere (mol/cm2)
        m_i : float
            Mass of the impactor (g)

        Returns
        -------
        N_H2O_steam : float
            moles/cm2 steam in the atmosphere.
        """

        T_0 = 300 # inital temperature
        C_H2O = self.partial_cp_cgs('H2O',T_0)
        C_w = C_H2O
        C_CO2 = self.partial_cp_cgs('CO2',T_0)
        C_N2 = self.partial_cp_cgs('N2',T_0)
        C_CO = self.partial_cp_cgs('CO',T_0)
        Q_w = 2.5e10 # latent heat of water vaporization ergs/g

        # convert to grams from mol/cm2
        M_H2O = N_H2O_ocean*self.gas.molecular_weights[self.gas.species_names.index('H2O')]*self.area
        M_CO2 = N_CO2*self.gas.molecular_weights[self.gas.species_names.index('CO2')]*self.area
        M_N2 = N_N2*self.gas.molecular_weights[self.gas.species_names.index('N2')]*self.area

        # heat capacity of dry atmosphere (erg/K)
        M_AC_A = M_CO2*C_CO2 + M_N2*C_N2

        # energy of impactor (ergs)
        E_i = 0.5*m_i*self.v_i**2

        # assume that the heated atm/steam is heated to T_prime
        M_steam = (E_i*self.impactor_energy_frac-M_AC_A*(self.T_prime-T_0))/(Q_w + (self.T_prime-T_0)*C_w)

        if M_steam > M_H2O: # Can't be more steam than water
            M_steam = M_H2O

        # convert to moles/cm2 of steam
        N_H2O_steam = M_steam/self.gas.molecular_weights[self.gas.species_names.index('H2O')]/self.area

        return N_H2O_steam

    def react_iron(self,N_H2O_steam,N_CO2,M_i):
        """Reacts iron from an impactor with oxidized gases in an atmosphere.
        Each mol of Fe sequesters one mol of O.

        Parameters
        ----------
        N_H2O_steam : float
            Steam in the atmosphere (mol/cm2)
        N_CO2 : float
            CO2 in the atmosphere (mol/cm2)
        M_i : float
            Mass of the impactor (g)

        Returns
        -------
        N_H2 : float
            H2 in the atmosphere (mol/cm2)
        N_H2O : float
            H2O in the atmosphere (mol/cm2)
        N_CO : float
            CO in the atmosphere (mol/cm2)
        N_CO2 : float
            CO2 in the atmosphere (mol/cm2)
        """
        Moles_Fe = self.Fe_react_frac*self.impactor_Fe_frac*M_i/56.
        N_Fe     = Moles_Fe/self.area
        xxx = N_Fe/(N_H2O_steam +N_CO2) #moles Fe/moles O
        if xxx<1:
            N_H2 = xxx*N_H2O_steam
            N_H2O = (1-xxx)*N_H2O_steam
            N_CO = xxx*N_CO2
            N_CO2 = (1-xxx)*N_CO2
        else:
            raise Exception('More Fe than O!')
        return N_H2, N_H2O, N_CO, N_CO2

class ImpactSolution():
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
            Tsurf, N, mubar, Psurf, Ntot, mix = prep_atm_first(y, self.stm)
        elif integ_type == 2:
            ind_H2O = self.stm.gas.species_names.index('H2O')
            Tsurf, N, P_H2O, mubar, Psurf, Ntot, mix = prep_atm_second(y, self.stm)
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

# rhs function for the first integration, when no condensation happens
def prep_atm_first(y, self):
    T = y[0]
    N = y[1:self.ngas+1]
    Ntot = np.sum(N)
    mix = N/Ntot
    self.gas.X = mix
    mubar = self.gas.mean_molecular_weight
    Psurf = Ntot*mubar*self.grav
    return T, N, mubar, Psurf, Ntot, mix

def prep_atm_first_catalyst(y, self):
    N_cat = y[self.ngas+1:]
    X_cat = N_cat/np.sum(N_cat)
    return X_cat

def rhs_first(t, y, dy, self):

    try:
        T, N, mubar, Psurf, Ntot, mix = prep_atm_first(y, self)
        self.gas.TPX = T, Psurf/10.0, mix
        if self.surface_catalyst:
            X_cat = prep_atm_first_catalyst(y, self)
            self.surf_phase.TPX = T,Psurf/10.0,X_cat
    except CanteraError:
        return 1
    
    # scale height
    Ha = const.N_avo*const.k_boltz*T/(mubar*self.grav)
    
    if self.surface_catalyst:
        all_rates = self.surf_phase.net_production_rates
        surf_rates = all_rates[:self.nsurf]*1.0e3/1.0e4 # mol/(cm2 metal * s)
        gas_rates = all_rates[self.nsurf:]*1.0e3/1.0e6 # mol/(cm3 * s)

        # surface rates
        dy[self.ngas+1:] = surf_rates*self.Ni_area*1.0e12 # pico mol/(cm2 Earth * s)
    else:
        gas_rates = self.gas.net_production_rates*1.0e3/1.0e6 # mol/(cm3 * s)
        
    # gas rates
    dy[1:self.ngas+1] = gas_rates*Ha # moles/cm2/s

    # climate
    Fir = OLR(T)
    dT_dt = -self.grav/(const.cp_H2O*Psurf)*Fir
    dy[0] = dT_dt
    
    return 0

# root finding function for first integration, where
# nothing is condensing
def stop_first(t, y, gout, self):
    try:
        T, N, mubar, Psurf, Ntot, mix = prep_atm_first(y, self)
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

def prep_atm_second(y, self):
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

def rhs_second(t, y, dy, self):

    try:
        T, N, P_H2O, mubar, Psurf, Ntot, mix = prep_atm_second(y, self)
        self.gas.TPX = T, Psurf/10.0, mix
    except CanteraError:
        return 1

    # scale height
    Ha = const.N_avo*const.k_boltz*T/(mubar*self.grav)
    
    # moles/cm3/s
    rx_rates = self.gas.net_production_rates*1.0e3/1.0e6
    # moles/cm2/s
    dN_dt = rx_rates*Ha
    
    # climate
    Fir = OLR(T)
    dT_dt = (-self.grav/(const.cp_H2O*Psurf))*Fir \
            *(1 + const.L_H2O**2*const.mu_H2O*P_H2O/(const.cp_H2O*Psurf*const.R*T**2))**-1

    dy[0] = dT_dt
    dy[1:self.ind_H2O+1] = dN_dt[:self.ind_H2O]
    dy[self.ind_H2O+1:] = dN_dt[self.ind_H2O+1:]
    
    return 0

def stop_second(t, y, gout, self):
    gout[0] = y[0] - self.T_stop
    return 0
    

