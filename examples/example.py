import numpy as np
from matplotlib import pyplot as plt
from EvolveAtm import Atmosphere as atm

# initial conditions (molecules/cm2)
Ninit_dict = {'H2':1.5e+26,
              'CO':9.2e+22,
              'CO2':6.3e+25,
              'CH4':6.5e+24,
              'N2':2.1e+25,
              'NH3':0.0,
              'HCN':0.0,
              'C2Hn':0.0,
              'Haze':0.0}

# integrate for 10 million years
tspan = [0,10e6*(60*60*24*365)]

out = atm.integrate(tspan,Ninit_dict) 

# plot results
plt.rcParams.update({'font.size': 20})
fig,ax = plt.subplots(1,1,figsize=[9,6])
yr = (365*24*60*60)
spec = ['H2','CO','CH4','CO2','N2','NH3']
for sp in spec:
    ax.plot(out['time']/yr/1e6,out[sp]*out['Ntot']/6.022e23,'-',label=sp)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(.1,1000)
ax.set_xlim(.03,1e1)
ax.set_ylabel('Column abundance (moles/cm$^2$)')
ax.set_xlabel('Time after impact (Myr)')
ax.legend(ncol=2)

ax2 = ax.twinx()
ax2.set_ylabel('HCN production\n(molecules cm$\mathrm{^{-2}}$ s$\mathrm{^{-1}}$)')  # we already handled the x-label with ax1
ax2.plot(out['time']/yr/1e6, out['dNHCN_dt'],'C7:',lw=4,label='HCN')
ax2.set_yscale('log')
ax2.legend()
ax2.set_ylim(1e7,2e10)
plt.show()
