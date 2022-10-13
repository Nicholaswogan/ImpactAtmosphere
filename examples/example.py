import numpy as np
from matplotlib import pyplot as plt
from ImpactAtmosphere_1 import EvolveAtm as atm

# initial conditions (moles/cm2)
Ninit_dict = {'H2': 298.9,
             'CO': 0.0562,
             'CO2': 22.9,
             'CH4': 0.0166,
             'N2': 35.96,
             'NH3': 0.0}

out = atm.integrate(Ninit_dict)

# plot results
plt.rcParams.update({'font.size': 20})
fig,ax = plt.subplots(1,1,figsize=[9,6])
yr = 365*24*60*60
spec = ['H2','CO','CH4','CO2','N2','NH3']
for sp in spec:
    ax.plot(out['time']/yr,out[sp]*out['Ntot'],'-',label=sp)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e-4,1000)
ax.set_xlim(1e2,1e6)
ax.set_ylabel('Column abundance (moles/cm$^2$)')
ax.set_xlabel('Time after impact (Myr)')
ax.legend(ncol=2)

ax2 = ax.twinx()
ax2.set_ylabel('HCN production\n(molecules cm$\mathrm{^{-2}}$ s$\mathrm{^{-1}}$)')
ax2.plot(out['time']/yr, out['dNHCN_dt'],'C7:',lw=4,label='HCN')
ax2.set_yscale('log')
ax2.legend()
ax2.set_ylim(1e7,2e10)
plt.show()
