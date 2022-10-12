import numpy as np
import cantera as ct
import os

from . import constants as const

class SteamAtm():

    def __init__(self,gas, T_prime = 2000, impactor_energy_frac = 0.5, \
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

        self.ind_H2O = self.gas.species_names.index('H2O')

        # Parameters
        self.T_prime = T_prime # Initial atmosphere temperature (K)
        self.impactor_energy_frac = impactor_energy_frac # fraction of kinetic energy used to heat the ocean and atmosphere
        self.Fe_react_frac = Fe_react_frac # fraction of the iron that reacts with atmosphere.
        self.impactor_Fe_frac = impactor_Fe_frac # Fe mass fraction of the impactor.
        self.v_i = v_i # velocity of impactor (cm/s)
        if self.surface_catalyst:
            self.Ni_area = Ni_area # cm^2 Ni surface / cm^2 column of the atmosphere

        # Settings
        self.rtol = 1e-12 # relative tolerance for CVODES integrator
        self.atol = 1e-21 # absolute tolerance for CVODES integrator
        self.dTemp = 10 # Temperature descritization (K)
        self.verbose = True

        # Constants
        self.grav = 982 # gravity of Earth (cm/s2)
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
        return sol.to_dict()

    def cooling_steam_atmosphere_1(self, N_H2O_ocean, N_init, P_init, X):
        """Integrates atmosphere from a very hot state (2000 K) down to where 
        H2O begins to condense.
        """
        
        # create solution object and append initial conditions
        sol = ImpactSolution(self.gas, self.grav)
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
            mubar_prev, Psurf_prev, Ntot_prev, Mtot_prev, mix_prev = prep_atmosphere(Tsurf_prev, N_prev, self.gas, self.grav)
            dtime = self.dtime_dTemp(Tsurf_prev, N_prev) # get the size of the time-step
            self.gas.TPX = Tsurf_prev,Psurf_prev/10.0,mix_prev
            
            # advance chemistry forward dtime using cantera
            r = ct.ConstPressureReactor(self.gas,energy='off')
            if self.surface_catalyst:
                self.surf_phase.TPX = Tsurf_prev,Psurf_prev/10.0,{'Ni':1}
                # advance material surface for 10 seconds while keeping
                # gas phase fixed. This makes the material properties
                # realistic and consistent with the gas
                self.surf_phase.advance_coverages(10,max_steps=1e7,rtol=self.rtol,atol=self.atol)
                # volume of a 1 cm2 column of the atmosphere. I multiply scale height
                # by 1.0 cm2
                Hscale_prev = const.N_avo*const.k_boltz*Tsurf_prev/(mubar_prev*self.grav)
                volume_column = Hscale_prev*1.0
                r.volume = volume_column/1e6 # convert to m3
                Ni_area_SI = self.Ni_area/1e2 # conver to m2
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
            # Only save if H2O is not supersaturated
            if Tsurf < const.T_crit_H2O:
                mubar, Psurf, Ntot, Mtot, mix = prep_atmosphere(Tsurf, N, self.gas, self.grav)
                P_H2O_sat = sat_pressure_H2O(Tsurf)
                P_H2O = Psurf*mix[self.ind_H2O]
                if P_H2O_sat > P_H2O:
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

        # tn = sol.time[-1]
        # Tsurf_prev = self.T_prime
        # Tsurf = self.T_prime
        # N_prev = N_init.copy()
        # N = N_init.copy()



        return None

    def dtime_dTemp(self, T, N):  
        mubar, Psurf, Ntot, Mtot, mix = prep_atmosphere(T, N, self.gas, self.grav)
        Fir = OLR(T)
        dT_dt = -self.grav/(const.cp_H2O*Psurf)*Fir
        dtime = -(1/dT_dt)*self.dTemp
        return dtime

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

    def steam_from_impact(self,N_H2O_ocean, N_CO2, N_N2, m_i):
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
    def __init__(self, gas, grav):
        self.gas = gas
        self.grav = grav
        self.time = np.empty((0,),np.float64)
        self.Tsurf = np.empty((0,),np.float64)
        self.mubar = np.empty((0,),np.float64)
        self.Psurf = np.empty((0,),np.float64)
        self.Ntot = np.empty((0,),np.float64)
        self.mix = {}
        for i,name in enumerate(self.gas.species_names):
            self.mix[name] = np.empty((0,),np.float64)

    def append(self, t, Tsurf, N):
        mubar, Psurf, Ntot, Mtot, mix = prep_atmosphere(Tsurf, N, self.gas, self.grav)
        # append
        self.time = np.append(self.time, t)
        self.Tsurf = np.append(self.Tsurf, Tsurf)
        self.mubar = np.append(self.mubar, mubar)
        self.Psurf = np.append(self.Psurf, Psurf)
        self.Ntot = np.append(self.Ntot, Ntot)
        for i,name in enumerate(self.gas.species_names):
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

def sat_pressure_H2O(T):
    return 1.0e6*np.exp(const.L_H2O*const.mu_H2O/const.R*(1/373.0 - 1/T)) # dynes/cm^2

# OLR for a steam atmosphere
def OLR(T):
    return 1.6e5 + 500.0*np.maximum(T-1200.0,0.0) # ergs/cm2/s

def prep_atmosphere(Tsurf, N, gas, grav):
    Ntot = np.sum(N)
    mix = N/Ntot
    gas.TPX = Tsurf, 1.0e5, mix
    mubar = gas.mean_molecular_weight
    Mtot = Ntot*mubar
    Psurf = Mtot*grav
    return mubar, Psurf, Ntot, Mtot, mix