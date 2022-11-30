import numpy as np
import cantera as ct
import os

class SteamAtmBase():
    """Base class for the two approaches to solving the steam atmosphere
    """
    
    def __init__(self, gas, T_prime, impactor_energy_frac, \
                 Fe_react_frac, impactor_Fe_frac, v_i, 
                 Ni_area):
        
        self.surface_catalyst = False
        if type(gas) == str:
            zahnle_path = os.path.dirname(os.path.realpath(__file__))+'/data/'
            ct.add_directory(zahnle_path)
            if "Methanation" in gas:
                # we have a surface catalyst
                self.surface_catalyst = True
                self.surf_phase = ct.Interface(gas,'Ni-surface')
                self.gas = self.surf_phase.adjacent["gas"]
            else:
                self.gas = ct.Solution(gas)
        else:
            self.gas = gas # cantera gas object
        self.gas.basis = 'mass'

        self.ngas = self.gas.n_total_species
        self.ngas_1 = 1
        if self.surface_catalyst:
            self.nsurf = self.surf_phase.n_total_species - self.ngas
        else:
            self.nsurf = 0
            self.nsurf_1 = 1+self.ngas

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

        # Planet
        self.grav = 981.0 # gravity of Earth (cm/s2)
        self.area = 5.1e18 # Area of Earth (cm2)

    ##########################################
    ### Stuff that happens right at impact ###
    ##########################################

    def initial_conditions(self, N_H2O_ocean, N_CO2, N_N2, M_i, N_CO_init, N_H2_init, N_CH4_init):
        """Determines the initial conditions for integration
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

        Ntot = np.sum(N_init)
        mix = N_init/Ntot
        self.gas.X = mix
        mubar = self.gas.mean_molecular_weight
        P_init = Ntot*mubar*self.grav
        
        self.gas.TPX = self.T_prime,P_init/10.0,mix
        self.gas.equilibrate('TP') # equilibrate

        # here we assume that Ntot has not changed much during equilibration
        N_init = self.gas.X*Ntot
        return N_init, P_init, self.gas.X

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
        N_Fe = Moles_Fe/self.area
        xxx = N_Fe/(N_H2O_steam +N_CO2) #moles Fe/moles O
        if xxx<1:
            N_H2 = xxx*N_H2O_steam
            N_H2O = (1-xxx)*N_H2O_steam
            N_CO = xxx*N_CO2
            N_CO2 = (1-xxx)*N_CO2
        else:
            raise Exception('More Fe than O!')
        return N_H2, N_H2O, N_CO, N_CO2

    ###################################
    ### Additional utility routines ###
    ###################################

    def partial_cp_cgs(self,spec,T):
        '''Finds specific heat capacity of spec in
        ergs/g/K (cgs units).
        '''
        self.gas.TP = T,1e5 # pressure doesn't matter
        index = self.gas.species_names.index(spec)
        cp = self.gas.partial_molar_cp[index]/self.gas.molecular_weights[index]*1e4
        return cp

    def make_equilibrium_solution(self, sol):
        out = {'time' : sol['time']}
        for key in self.gas.species_names:
            out[key] = np.empty(len(sol['time']))

        for i in range(len(sol['time'])):
            T = sol['Tsurf'][i]
            P = sol['Psurf'][i]/10.0 # Pa
            X = np.empty(self.gas.n_total_species)
            for j,key in enumerate(self.gas.species_names):
                X[j] = sol[key][i]
            self.gas.TPX = T,P,X
            self.gas.equilibrate("TP")
            
            for j,key in enumerate(self.gas.species_names):
                out[key][i] = self.gas.X[j]
        return out

    def init_for_integrate(self, sol):
        Ninit_dict = {}
        species = ['H2','CO','CO2','CH4','N2','NH3']
        for spec in species:
            if spec in sol:
                Ninit_dict[spec] = sol[spec][-1]*sol['Ntot'][-1]
            else:
                Ninit_dict[spec] = 0.0
        Ninit_dict['NH3'] = 0.0
        return Ninit_dict
        
    def dry_end_atmos(self, sol_stm, clip = True):
        Ntot = sol_stm['Ntot'][-1]
        Ntot_dry = Ntot - Ntot*sol_stm['H2O'][-1]

        sol_dry = {}
        sol_dry['Ntot'] = Ntot_dry

        mubar_dry = 0.0
        for i,sp in enumerate(self.gas.species_names):
            if sp =='H2O':
                sol_dry['H2O'] = 0.0
            else:
                if clip:
                    mix = np.maximum(sol_stm[sp][-1],0.0)
                else:
                    mix = sol_stm[sp][-1]
                sol_dry[sp] = mix*Ntot/Ntot_dry
                mubar_dry += self.gas.molecular_weights[i]*sol_dry[sp]
        Psurf_dry = Ntot_dry*mubar_dry*self.grav
        sol_dry['Psurf'] = Psurf_dry
        sol_dry['mubar'] = mubar_dry
        
        return sol_dry

