import numpy as np
from photochem import photochem

def evolve(N_in,fH2O):
    LH2O,LH2,LCO,LCO2,LCH4,LN2,LNH3,LHCN,LC2Hn,LHaze = \
    0,1,2,3,4,5,6,7,8,9
    N_avo = 6.02e23
    names = ['H2O','H2','CO','CO2','CH4',\
             'N2','NH3','HCN','C2Hn','Haze']

    # N_in is moles/cm2
    # N_in1 is molecules/cm2
    N_in1 = np.zeros(10)
    for i,name in enumerate(names):
        N_in1[i] = N_in[name]*N_avo

    Time,T_surf,pressure,mu,N_tot,f_tot,dNdt = photochem(N_in1,fH2O)
    ind = np.count_nonzero(T_surf)-1 # where the simulation stopped
    Time,T_surf,pressure,mu,N_tot,f_tot,dNdt = \
    Time[:ind],T_surf[:ind],pressure[:ind],mu[:ind],\
    N_tot[:ind]/N_avo,f_tot[:ind],dNdt[:ind]


    # convert output to dict
    output = {}
    output['time'] = Time
    output['temp'] = T_surf
    output['Pres'] = pressure
    output['mu'] = mu
    output['Ntot'] = N_tot
    # mixing ratios
    for i,name in enumerate(names):
        output[name] = f_tot[:,i]

    # production rates
    for i,name in enumerate(names):
        output['dN_'+name+'dt'] = dNdt[:,i]

    return output
