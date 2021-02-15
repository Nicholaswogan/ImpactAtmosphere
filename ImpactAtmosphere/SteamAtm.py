import numpy as np
import cantera as ct
import sys
import os

from .EvolveAtm import integrate

class SteamAtm():

    def __init__(self,gas,T_prime = 2000, impactor_energy_frac = 0.5, \
        Fe_react_frac = 1, impactor_Fe_frac = 0.33, v_i = 17e5 ):

        if type(gas) == str:
            zahnle_path = os.path.dirname(os.path.realpath(__file__))+'/data/'
            ct.add_directory(zahnle_path)
            self.gas = ct.Solution(gas)
        else:
            self.gas = gas # cantera gas object
        self.gas.basis = 'mass'

        # Parameters
        self.T_prime = T_prime # Initial atmosphere temperature (K)
        self.impactor_energy_frac = impactor_energy_frac # fraction of kinetic energy used to heat the ocean and atmosphere
        self.Fe_react_frac = Fe_react_frac # fraction of the iron that reacts with atmosphere.
        self.impactor_Fe_frac = impactor_Fe_frac # Fe mass fraction of the impactor.
        self.v_i = v_i # velocity of impactor (cm/s)

        # Settings
        self.rtol = 1e-12 # relative tolerance for CVODES integrator
        self.atol = 1e-21 # absolute tolerance for CVODES integrator
        self.dTemp = 10 # Temperature descritization (K)
        self.verbose = True

        # Constants
        self.grav = 982 # gravity of Earth (cm/s2)
        self.area = 5.1e18 # Area of Earth (cm2)
        self.k_B = ct.boltzmann*1e7 # boltzman constant (cgs units)
        self.Navo = ct.avogadro*1e-3 # avagadros's number (molecules/mol)

    def impact(self,N_H2O_ocean,N_CO2,N_N2,M_i):
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
        self.initial_conditions(N_H2O_ocean,N_CO2,N_N2,M_i)

        self.gas.TPX = self.T_prime,P_init*0.1,X

        # Keep track of these things
        Temperature = np.array([self.T_prime]) # K
        Column = np.array([np.sum(N_init)]) # mol/cm2
        mubar = np.array([self.gas.mean_molecular_weight]) # g/mol
        seconds = np.array([0]) # s
        Mix = np.array([X],)
        Mass = Column[0]*self.gas.mean_molecular_weight # g/cm2

        while Temperature[-1]>300:
            if self.verbose:
                print('Atmospheric Temperature = '+"{:<8}".format('%.1f'%Temperature[-1]+' K'),end='\r')
            # compute the timestep
            dtime, break_flag = self.dtime_dTemp(Temperature[-1], \
                                P_init, \
                                mubar[-1], \
                                Mix[-1,self.gas.species_names.index('H2O')]*P_init)

            if break_flag:
                break

            # advance to dtime using cantera
            self.gas.TPX = Temperature[-1],P_init*.1,Mix[-1]
            r = ct.ConstPressureReactor(self.gas,energy='off') # no T change during reaction
            reactorNetwork = ct.ReactorNet([r])
            reactorNetwork.rtol = self.rtol
            reactorNetwork.atol = self.atol
            t = 0.0
            while(t < dtime):
                t = reactorNetwork.step()
            reactorNetwork.advance(dtime)
            self.gas.X = self.gas.X.clip(min=0) # force positive

            # save results
            Temperature = np.append(Temperature,Temperature[-1]-self.dTemp)
            mubar = np.append(mubar,self.gas.mean_molecular_weight)
            Column = np.append(Column,Mass/mubar[-1])
            seconds = np.append(seconds,seconds[-1]+dtime)
            Mix = np.append(Mix,[self.gas.X],axis=0)

        if not break_flag:
            sys.exit('Integration failed. Water never started to condense')
        if self.verbose:
            print("{:<40}".format('Integration successful.'))

        # build dictionary of output
        solution = {'Tsurf': Temperature,\
                    'Psurf': P_init/1e6,\
                    'Ntot': Column,\
                    'mubar' : mubar,\
                    'time' : seconds}
        for i,name in enumerate(self.gas.species_names):
            solution[name] = Mix[:,i]
        return solution

    def init_for_integrate(self,solution):
        Ninit_dict = {}
        species = ['H2','CO','CO2','CH4','N2','NH3']
        for spec in species:
            Ninit_dict[spec] = solution[spec][-1]*solution['Ntot'][-1]
        Ninit_dict['NH3'] = 0.0
        return Ninit_dict

    def impact_integrate(self,N_H2O_ocean,N_CO2,N_N2,M_i,**kwargs):
        solution = self.impact(N_H2O_ocean,N_CO2,N_N2,M_i)
        N_init_dict = self.init_for_integrate(solution)
        out = integrate(Ninit_dict,**kwargs)
        return out

    def initial_conditions(self,N_H2O_ocean,N_CO2,N_N2,M_i):
        """Determines the initial conditions for method `impact`.
        """

        # Determine how much steam is made
        N_H2O_steam = self.steam_from_impact(N_H2O_ocean,N_CO2,N_N2,M_i)
        # React the atmosphere with impactor Fe
        N_H2,N_H2O,N_CO,N_CO2 = self.react_iron(N_H2O_steam,N_CO2,M_i)
        # set as initial conditions
        N_init = np.zeros(len(self.gas.species_names))
        N_init[self.gas.species_names.index('H2')] = N_H2
        N_init[self.gas.species_names.index('H2O')] = N_H2O
        N_init[self.gas.species_names.index('CO')] = N_CO
        N_init[self.gas.species_names.index('CO2')] = N_CO2
        N_init[self.gas.species_names.index('N2')] = N_N2

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

    def dtime_dTemp(self,T_in, p, mu, pH2O):
        """Determines the amount of time (s) elapsing between temperature steps
        (self.dTemp). Also determines if steam is beginning to condense.

        Parameters
        ----------
        T_in : float
            Temperature (K)
        p : float
            Pressure (dynes)
        mu : float
            mean molecular weight of atmosphere (g/mol)
        pH2O : float
            Partial pressure of H2O (dynes)

        Returns
        -------
        dtime : float
            Time in seconds it takes to atmosphere to cool from T_in to T_in-self.dTemp
        break_flag : bool
            If True, then water is beginning to condense.
        """
        gamma = 1.33 # ???
        L_v = 2.5e10 # latent heat of water vaporization erg/K
        AA,TB = 3.88e11,4800

        T     = T_in - 0.5*self.dTemp
        T_out = T_in - self.dTemp

        Hscale = self.Navo*self.k_B*T/(mu*self.grav)   # scale height with fixed mu
        F_ir = 1.6e5 + 500.*np.max([T-1200.,0.])
        bracket = gamma/(gamma-1)*p*Hscale/T
        dtime = bracket/F_ir*self.dTemp

        # determine if water is beginning to saturate
        psat = AA*np.exp(-TB/T_out)
        if psat < pH2O:
            break_flag = True
        else:
            break_flag = False

        return dtime, break_flag

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
            sys.exit('More Fe than O!')
        return N_H2, N_H2O, N_CO, N_CO2

    def diameter(self,mm): # in g
        m = mm/1e3 # to kg
        rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
        return 2*((3/(4*np.pi))*m/rho)**(1/3) # km

    def mass(self,D): #km
        rho = 3.0e12 # kg/km3 (Morbidelli et al. 2012)
        return rho*(4/3)*np.pi*(D/2)**3*1e3 # mass in g
