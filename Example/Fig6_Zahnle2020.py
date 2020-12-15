import numpy as np
from matplotlib import pyplot as plt
from EvolveAtm.photochemistry import EvolveAtm

ev = EvolveAtm()
t = 0 # start at t=0

# Atmosphere in molecules/cm2
N = np.array([1.68675180e+28, 1.53771633e+26, 9.21714507e+22, 6.30403649e+25,
       6.55363284e+24, 2.18545381e+25, 9.36588326e+22, 1.00000000e+10,
       1.00000000e+10, 1.00000000e+10])

t = [0,10e6*(365*24*60*60)]
t_eval = np.linspace(0,10e6*(365*24*60*60),1000)
output = ev.integrate(t,N,out_dict=True,t_eval=t_eval)

t = output['time']
H2 = output['H2']*output['press']
CO = output['CO']*output['press']
CH4 = output['CH4']*output['press']
CO2 = output['CO2']*output['press']
H2O = output['H2O']*output['press']
NH3 = output['NH3']*output['press']
N2 = output['N2']*output['press']

Haze=output['dNHaze_dt']
HCN=output['dNHCN_dt']


yr = (365*24*60*60)
plt.rcParams.update({'font.size': 15})
fig,ax = plt.subplots(1,1,figsize=[9,6])
ax.plot(t/yr/1e6,H2,'-',label='H2')
ax.plot(t/yr/1e6,CO,label='CO')
ax.plot(t/yr/1e6,CH4,label='CH4')
ax.plot(t/yr/1e6,CO2,label='CO2')
ax.plot(t/yr/1e6,N2,'C5',label='N2')
ax.plot(t/yr/1e6,NH3,'C6',label='NH3')


ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e-2,10)
ax.set_xlim(0.03,2)
ax.set_ylabel('Partial pressure (bar)')
ax.set_xlabel('Time after quench (Myr)')
ax.legend(ncol=2)

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('HCN, Haze Production\n(molecules cm$\mathrm{^{-2}}$ s$\mathrm{^{-1}}$)')  # we already handled the x-label with ax1
ax2.plot(t/yr/1e6, HCN,'C4:',lw=4,label='HCN')
ax2.plot(t/yr/1e6, Haze,'k:',lw=4,label='Haze')
ax2.set_yscale('log')
ax2.set_ylim(1e7,2e10)
ax2.legend()

plt.savefig("Fig6_Zahnle2020.pdf",bbox_inches='tight')
plt.show()
