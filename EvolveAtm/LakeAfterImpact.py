from EvolveAtm import *
# this is an extension of EvolveAtm


class LakeAfterImpact:

    def __init__(self):
        self.species = ['H2O','H2','CO','CO2','CH4','N2','NH3', \
                        'HCN', 'C2Hn', 'Haze','HCONH2','HCOOH']

    def parameters(self,T_surf,pH_ocean,T_lake,pH_lake,Kzz = 1.e5,\
                   HCNalt = 60.0e5, nz = 60, T_trop = 180., \
                   P_trop = 0.1, vo = 1.2e-5, zs = 100.e2, zd = 4000.e2, \
                   fH2O_trop = 1.0e-6):

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
        HCONH2, HCOOH = y[-2:]

        # the right-hand-side of the atmosphere
        dNdt, pressure, T_surf, tau_uv, N_H2O, mubar = atmos.rhs_verbose(t,N)

        # compute diffusion of HCN to the surface
        alt, fHCN = diffuse(dNdt[self.species.index('HCN')], Ts = self.T_surf, \
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
        dHCONH2_dt = self.ktot_lake*sumHCN - self.kform_lake*HCONH2
        dHCOOH_dt = self.kform_lake*HCONH2

        # right-hand-side of ODEs
        RHS = np.append(dNdt,np.array([dHCONH2_dt,dHCOOH_dt]))

        return RHS, pressure, T_surf, tau_uv, N_H2O, mubar, PHCN, HCN_aq, CN, sumHCN

    def rhs(self,t,y):
        RHS, pressure, T_surf, tau_uv, N_H2O, mubar, PHCN, HCN_aq, CN, sumHCN = self.rhs_verbose(t,y)
        return RHS

    def integrate(self,tspan,init_dict,out_dict=True,method = "LSODA",**kwargs):
        init = np.zeros(len(self.species))
        init[0] = 1. # H2O doesn't matter
        init[1:] = np.array([init_dict[key] for key in self.species[1:]])

        # set tau uv to 100. then integrate
        atmos.tau_uv_init = 100.
        sol = solve_ivp(self.rhs,tspan,init,method = method,**kwargs)

        # check if solver succeeded
        if not sol.success:
            sys.exit('ODE integration failed')

        # only take positive solution (it should be almost positive)
        yvals = np.clip(sol.y.T,0,np.inf)
        tvals = sol.t

        if not out_dict:
            return tvals,yvals
        elif out_dict:
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
            out = {}
            out['Ntot'] = np.sum(yvals[:,:-5],axis=1)-\
                          yvals[:,self.species.index('H2O')] + N_H2O
            out['press'] = pressure
            out['Tsurf'] = T_surf
            out['tau_uv'] = tau_uv
            out['mubar'] = mubar
            out['time'] = tvals
            out['dNHCN_dt'] = dydt[:,self.species.index('HCN')]
            out['dNHaze_dt'] = dydt[:,self.species.index('Haze')]
            out['H2O'] = N_H2O/out['Ntot']
            out['NH2O_strat'] = yvals[:,self.species.index('H2O')]
            out['PHCN'] = PHCN
            out['HCN_aq_lake'] = HCN_aq
            out['CN_lake'] = CN
            out['sumHCN_lake'] = sumHCN
            out['HCONH2_lake'] = yvals[:,-2]
            out['HCOOH_lake'] = yvals[:,-1]
            for i in range(1,len(self.species)-2):
                    out[self.species[i]] = yvals[:,i]/out['Ntot']

            return out
