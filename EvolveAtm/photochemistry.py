import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
import sys
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
import os


class EvolveAtm:
    # constants
    area = 5.1e18 #cm2
    g = 982.0 #cm/s2

    LH2O,LH2,LCO,LCO2,LCH4,LN2,LNH3, LHCN, LC2Hn, LHaze = \
    0,1,2,3,4,5,6,7,8,9
    species = ['H2O','H2','CO','CO2','CH4','N2','NH3', 'HCN', 'C2Hn', 'Haze']
    # masses in g/mol
    mass = 1.0e3*np.array([0.016,0.002,0.028,0.044,0.016,0.028,0.017,0.027,0.026,0.099])
    N_avo = 6.022e23
    Nsun = 7
    Nsp = 7
    Nspecies = 10
    A_escape, B_escape = [0.723*2e12, 0.006]
    Wolf_Toon = 3e10
    eta = 0.14
    p_trop = 0.1 # bar
    T_trop = 180
    fH2O = 1.0e-6
    Sun = np.array([30,20,10,10,10, 4,2])
    SSX = 1.0
    others3 = 1e6*1e4   #  scale height * density
    others4 = 1e6*1e4   #  scale height * density
    tau_uv_init = 100.
    data_path = os.path.dirname(os.path.realpath(__file__))+'/'
    # load other stuff
    k = np.loadtxt(data_path+'react_rates.txt')
    sigma = np.loadtxt(data_path+'sigma.txt')
    Flux0 = np.loadtxt(data_path+'Flux0.txt')
    Flux = Flux0*Sun*SSX

    def ode_eq_all(self,t,y,tau_uv):
        N = y

        ###### water vapor #######
        # first thing to do is to calculate how much water is in
        # the atmosphere... sorta correct some stuff
        N_dry = np.sum(N)-N[self.LH2O]
        #average g/mol in dry atmosphere
        mu_dry =(np.sum(N*self.mass)-N[self.LH2O]*self.mass[self.LH2O])/N_dry

        # This little routine self consistently calculates
        # surface pressure, temperature, and the H2O pressure
        # at the surface.
        p_surf = (mu_dry*N_dry*self.g/self.N_avo)
        T_s = self.T_trop*(p_surf/(self.p_trop*1e6))**(self.eta)
        p_H2O_surf  = 1e6*np.exp(5000./373.-5000./T_s )

        def system(y):
            T_s,P_s,P_H2O = y
            return(T_s-self.T_trop*(P_s/(self.p_trop*1e6))**self.eta,\
                   P_s-((N_dry + P_H2O*self.N_avo/(((P_s-P_H2O)*mu_dry+P_H2O*self.mass[self.LH2O])/P_s)/self.g)\
                        *(((P_s-P_H2O)*mu_dry+P_H2O*self.mass[self.LH2O])/P_s)*self.g/self.N_avo),\
                   P_H2O-1e6*np.exp(5000/373.-5000/T_s))

        sol = root(system,[T_s,p_surf,p_H2O_surf],method='lm',options={'maxiter':3000})
        error = np.linalg.norm(system(sol['x']))
        tol = 1e-7
        if sol.success==False:
            sys.exit('Failed to find surface pressure and temperature.')


        # name a bunch of stuff
        T_surf, pressure,p_H2O = sol['x']
        pressure = pressure/1e6
        #T_surf, pressure =  [299.69529487629842  ,      3.8148860339577828 ]
        p_H2O = np.exp(5000/373-5000/T_surf)
        p_dry = pressure-p_H2O
        mu = (p_dry*mu_dry + p_H2O*self.mass[self.LH2O])/pressure
        N_H2O = pressure*self.N_avo*1e6/mu/self.g - N_dry

        # water partial pressure in tropopause
        p_H2O_trop  = self.fH2O*self.p_trop

        # fraction of dry atmosphre above tropopause
        fraction = (self.p_trop-p_H2O_trop)/p_dry
        NX = N_dry*fraction
        #calculate statospheric H2O
        N_H2O_strat = p_H2O_trop/self.p_trop*NX

        # total stratospheric column??? Not sure of purpose of this
        N_strat  = NX + N_H2O_strat

        # statospheric H2O for photolysis purposes
        NH2Ox = p_H2O_trop/self.p_trop*N_dry
        N[self.LH2O] = NH2Ox

        ###### end water vapor #######

        ###### hydrogen escape
        N_bar = np.sum(N)-N[self.LH2]
        mubar = (np.sum(N[0:self.Nsp]*self.mass[0:self.Nsp])-N[self.LH2]*self.mass[self.LH2])/N_bar
        N_t = np.sum(N)

        #escape bit
        BBB = ((self.mass[self.LCO2]-self.mass[self.LH2])/(mubar-self.mass[self.LH2]))**2

        # units of molecules/s
        fH2 = N[self.LH2]/N_t
        H2_esc = -self.A_escape*self.Sun[1]*self.SSX*(fH2/(1+2.3*fH2))\
        /np.sqrt(1 + self.B_escape*(self.Sun[1]*self.SSX)**2*BBB)

        ##### Photolysis #####
        # this is the optical thickness from the haze
        tau = np.ones(self.Nsp)*tau_uv *N_t/NX

        # this is the optical thickness from molecules
        for L in range(0,self.Nsun):
            for j in range(0,self.Nsp):
                tau[L] = tau[L]+self.sigma[L,j]*N[j]

        # ion production
        Phi_ion = N[0:self.Nsp]*self.sigma[0,:]*self.Flux[0]/4/tau[0]


        Phi = np.zeros(self.Nsp)
        ### CO2 photolysis
        Phi[self.LCO2] = np.sum(self.Flux[1:6]/4*self.sigma[1:6,self.LCO2]*N[self.LCO2]/tau[1:6])
        Phi[self.LCO2] += 0.2*Phi_ion[self.LCO2]


        ### CH4 photolysis
        Phi[self.LCH4] = np.sum(self.Flux[0:4]/4*self.sigma[0:4,self.LCH4]*N[self.LCH4]/tau[0:4])

        #ion production
        sum_ions = np.sum(N[0:6])

        CH4_ion_breakup = Phi_ion[self.LH2]/sum_ions\
                        + Phi_ion[self.LN2]*(N[self.LCH4]/(sum_ions - N[self.LN2])\
                        + N[self.LCO]/(sum_ions -N[self.LN2])*N[self.LCH4]\
                        /(sum_ions -N[self.LN2] -N[self.LCO])\
                        + N[self.LCO]/(sum_ions -N[self.LN2])*N[self.LCO2]\
                        /(sum_ions -N[self.LN2] -N[self.LCO])\
                        * N[self.LCH4]/(N[self.LH2]+N[self.LCH4]+N[self.LH2O])\
                        + N[self.LH2O]/(sum_ions - N[self.LN2])\
                        * N[self.LCH4]/(N[self.LCH4]+N[self.LH2O])  )\
                        + Phi_ion[self.LCO]*( N[self.LCH4]/(N[self.LCO2]+N[self.LCH4]+N[self.LH2O])\
                        + N[self.LCO2]/(N[self.LCO2]+N[self.LCH4]+N[self.LH2O])\
                        * N[self.LCH4]/(N[self.LH2]+N[self.LCH4]+N[self.LH2O]) )\
                        + Phi_ion[self.LCO2]*N[self.LCH4]\
                        /(N[self.LH2]+N[self.LCH4]+N[self.LH2O])\
                        + Phi_ion[self.LH2O]*N[self.LCH4]/(N[self.LH2]+N[self.LCH4])

        Phi[self.LCH4] += CH4_ion_breakup

        ### H2O photolysis
        Phi[self.LH2O] = np.sum(self.Flux[1:6]/4*self.sigma[1:6,self.LH2O]*N[self.LH2O]/tau[1:6])

        ### N2 photolysis
        Phi[self.LN2] = np.sum(self.Flux[1:3]/4*self.sigma[1:3,self.LN2]*N[self.LN2]/tau[1:3])

        ### NH3 photolysis
        Phi[self.LNH3] = np.sum(self.Flux[0:7]/4*self.sigma[0:7,self.LNH3]*N[self.LNH3]/tau[0:7])

        ##### end Photolysis #####

        ##### Photochemistry #####
        ### O1D
        sumO1D = np.sum(self.k[1,0:6]*N[0:6])
        O1D_H2  = self.k[1,self.LH2]*N[self.LH2]/sumO1D
        O1D_CH4 = self.k[1,self.LCH4]*N[self.LCH4]/sumO1D

        ### O1D > OH
        O1D_OH = ( self.k[1,self.LH2]*N[self.LH2] + self.k[1,self.LCH4]*N[self.LCH4]\
                     + 2*self.k[1,self.LH2O]*N[self.LH2O] )/sumO1D

        ### OH
        sumOH = np.sum(self.k[2,0:5]*N[0:5])+self.k[2,5]*self.others3 # wtf
        OH_CH4 = self.k[2,self.LCH4]*N[self.LCH4]/sumOH
        OH_H2  = self.k[2,self.LH2]*N[self.LH2]/sumOH
        OH_CO  = self.k[2,self.LCO]*N[self.LCO]/sumOH

        ### N > HCN
        sumN2D = np.sum(self.k[0,0:6]*N[0:6])
        N_HCN = (self.k[0,self.LCH4]*N[self.LCH4] + self.k[0,self.LCO]*N[self.LCO]\
                    + self.k[0,self.LN2]*N[self.LN2] )/sumN2D

        ### O + CO > CO2
        sumO = np.sum(self.k[3,0:5]*N[0:5]) + self.k[3,5]*self.others4
        O_CH4 = self.k[3,self.LCH4]*N[self.LCH4]/sumO
        O_H2  = self.k[3,self.LH2]*N[self.LH2]/sumO
        O_CO  = self.k[3,self.LCO]*N[self.LCO]/sumO

        ##### end Photochemistry #####

        ##### Budgets #####
        # CH4 oxidiation
        dCH4dt_ox = -Phi[self.LCO2]*O1D_CH4 -Phi[self.LH2O]*OH_CH4

        # NH3 > HCN
        NH3_HCN =  (Phi[self.LCH4]-dCH4dt_ox)\
                    / (Phi[self.LCH4]-dCH4dt_ox + \
                       Phi[self.LH2O] +Phi[self.LCO2] )

        # CH4
        Phi_geo = np.zeros(self.Nsp)
        Phi_geo[self.LCH4] = 0.0   # outgassing of H2
        dNdt = np.zeros(self.Nspecies)
        dNdt[self.LCH4] = -Phi[self.LCH4] + dCH4dt_ox + Phi_geo[self.LCH4]

        # HCN
        dNdt[self.LHCN] = 2*Phi[self.LN2]*N_HCN*\
                        (Phi[self.LCH4]-dCH4dt_ox)\
                        /( Phi[self.LCH4]-dCH4dt_ox + 2*Phi[self.LN2]*N_HCN )\
                        + Phi[self.LNH3]*NH3_HCN

        # C2Hn
        dNdt[self.LC2Hn] = 0.5*(Phi[self.LCH4]-dCH4dt_ox)**3\
                        / (Phi[self.LCH4]-dCH4dt_ox + Phi[self.LH2O] \
                           +Phi[self.LCO2] )**2

        # Haze
        dNdt[self.LHaze] = (Phi[self.LCH4]-dCH4dt_ox) *\
                    ( (Phi[self.LCH4]-dCH4dt_ox)\
                    /(Phi[self.LCH4]-dCH4dt_ox + Phi[self.LH2O] +Phi[self.LCO2]) )**5

        # CO2
        Phi_geo[self.LCO2] = 0.0*N[self.LCO2]   # net subduction of CO2.
        dNdt[self.LCO2] = -Phi[self.LCO2] + Phi[self.LCO2]*O1D_OH*OH_CO\
                    + Phi[self.LH2O]*OH_CO + Phi[self.LCO2]*(1-O1D_OH)*O_CO\
                    - Phi_geo[self.LCO2]


        # CO
        dNdt[self.LCO] = Phi[self.LCO2]-(dCH4dt_ox-Phi[self.LCH4])\
                - dNdt[self.LHaze] - Phi[self.LCO2]*O1D_OH*OH_CO\
                - Phi[self.LH2O]*OH_CO - Phi[self.LCO2]*(1-O1D_OH)*O_CO\
                - dNdt[self.LHCN]

        # H2
        Phi_geo[self.LH2] = 0.0   # outgassing of H2

        dNdt[self.LH2] = H2_esc - 2*(dCH4dt_ox-Phi[self.LCH4])\
                + 2*dNdt[self.LCO2] + dNdt[self.LCO] - dNdt[self.LHaze] - dNdt[self.LHCN]\
                + Phi_geo[self.LH2] + 3*Phi[self.LNH3]*NH3_HCN


        # H2O
        dNdt[self.LH2O] = 0.0     # assumed constant

        # NH3
        dNdt[self.LNH3] = -Phi[self.LNH3]     # assumed lost

        # N2
        dNdt[self.LN2] = -0.5*dNdt[self.LHCN] +0.5*Phi[self.LNH3]*(1-NH3_HCN)

        return dNdt,pressure,T_surf,N_H2O

    def ode_eq(self,t,y,tau_uv):
        dNdt,pressure,T_surf,N_H2O = self.ode_eq_all(t,y,tau_uv)
        return dNdt

    def rhs(self,t,N):
        for i in range(0,len(N)):
            N[i] = np.max([0,N[i]])

        def func(tau_uv_log):
            dNdt = self.ode_eq(t,N,np.exp(tau_uv_log))
            return np.exp(tau_uv_log)-10.*((dNdt[self.LHaze]+dNdt[self.LHCN])/self.Wolf_Toon)**0.8
        sol = root(func,np.log(self.tau_uv_init))

        tau_uv = np.exp(sol.x[0])
        self.tau_uv_init = np.exp(sol.x[0])

        if tau_uv!=tau_uv:
            tau_uv = 0

        return self.ode_eq(t,N,tau_uv)

    def rhs_verbose(self,t,N):
        for i in range(0,len(N)):
            N[i] = np.max([0,N[i]])

        def func(tau_uv_log):
            dNdt = self.ode_eq(t,N,np.exp(tau_uv_log))
            return np.exp(tau_uv_log)-10.*((dNdt[self.LHaze]+dNdt[self.LHCN])/self.Wolf_Toon)**0.8
        sol = root(func,np.log(self.tau_uv_init))

        tau_uv = np.exp(sol.x[0])
        self.tau_uv_init = np.exp(sol.x[0])

        if tau_uv!=tau_uv:
            tau_uv = 0
        dNdt,pressure,T_surf,N_H2O = self.ode_eq_all(t,N,tau_uv)

        return dNdt,pressure,T_surf,tau_uv,N_H2O

    def integrate(self,tspan,Ninit,out_dict=False,t_eval = None,method='LSODA'):
        sol = solve_ivp(self.rhs,tspan,Ninit,\
                        method=method,t_eval = t_eval)

        # only take positive solution
        N_vals = np.clip(sol.y.T,1,np.inf)
        tvals = sol.t

        # check if solver succeeded
        if not sol.success:
            sys.exit('ODE integration failed')

        if not out_dict:
            return N_vals
        elif out_dict:
            dNdt = np.zeros(N_vals.shape)
            pressure = np.zeros(len(tvals))
            T_surf = np.zeros(len(tvals))
            tau_uv = np.zeros(len(tvals))
            N_H2O = np.zeros(len(tvals))
            self.tau_uv_init = 100
            for i in range(0,len(tvals)):
                dNdt[i],pressure[i],T_surf[i],tau_uv[i],N_H2O[i] = \
                self.rhs_verbose(tvals[i],N_vals[i])

            # output as dict
            keys = self.species+['press','Tsurf','tau_uv','Ntot']
            out = {}
            out['Ntot'] = np.sum(N_vals[:,:-3],axis=1)-\
                          N_vals[:,self.LH2O] + N_H2O
            out['press'] = pressure
            out['Tsurf'] = T_surf
            out['tau_uv'] = tau_uv
            out['time'] = tvals
            out['dNHCN_dt'] = dNdt[:,self.LHCN]
            out['dNHaze_dt'] = dNdt[:,self.LHaze]
            out['H2O'] = N_H2O/out['Ntot']
            out['NH2O_strat'] = N_vals[:,self.LH2O]

            for i in range(len(self.species)):
                if self.species[i] != 'H2O':
                    out[self.species[i]] = N_vals[:,i]/out['Ntot']

            return out

    def integrate_CVODE(self,t0,tfinal,Ninit,out_dict=False,ncp = 1000, ncp_list = None):
        mod = Explicit_Problem(self.rhs, Ninit, t0)
        sim = CVode(mod)
        sim.rtol = 1e-4
        sim.atol = 1e-6
        sim.verbosity = 50

        t,N_vals = sim.simulate(tfinal,ncp,ncp_list = ncp_list)


        N_vals = np.clip(N_vals,1,np.inf)
#         N_vals = sol.y.T
        tvals = np.array(t)

        if not out_dict:
            return N_vals
        elif out_dict:
            dNdt = np.zeros(N_vals.shape)
            pressure = np.zeros(len(tvals))
            T_surf = np.zeros(len(tvals))
            tau_uv = np.zeros(len(tvals))
            N_H2O = np.zeros(len(tvals))
            self.tau_uv_init = 100
            for i in range(0,len(tvals)):
                dNdt[i],pressure[i],T_surf[i],tau_uv[i],N_H2O[i] = \
                self.rhs_verbose(tvals[i],N_vals[i])

            # output as dict
            keys = self.species+['press','Tsurf','tau_uv','Ntot']
            out = {}
            out['Ntot'] = np.sum(N_vals[:,:-3],axis=1)-\
                          N_vals[:,self.LH2O] + N_H2O
            out['press'] = pressure
            out['Tsurf'] = T_surf
            out['tau_uv'] = tau_uv
            out['time'] = tvals
            out['dNHCN_dt'] = dNdt[:,self.LHCN]
            out['dNHaze_dt'] = dNdt[:,self.LHaze]
            out['H2O'] = N_H2O/out['Ntot']
            out['NH2O_strat'] = N_vals[:,self.LH2O]

            for i in range(len(self.species)):
                if self.species[i] != 'H2O':
                    out[self.species[i]] = N_vals[:,i]/out['Ntot']

            return out
