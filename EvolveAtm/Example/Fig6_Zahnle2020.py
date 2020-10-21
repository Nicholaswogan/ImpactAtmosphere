from EvolveAtm.EvolveAtm import evolve
import numpy as np
from matplotlib import pyplot as plt

yr = (365*24*60*60)

# Figure 6 case in Zahnle et al 2020
# mol/cm2
N_in = {}
N_in['CO'] = 0.15305787230451487
N_in['CO2'] = 104.68343563325219
N_in['H2'] = 255.34977197875762
N_in['H2O'] = 28009.827339782918
N_in['CH4'] = 10.882817733102817
N_in['N2'] = 36.291162496361459
N_in['NH3'] = 0.15552778586304855
N_in['HCN'] = 0.0
N_in['C2Hn'] = 0.0
N_in['Haze'] = 0.0


fH2O = 1.0e-6

output = evolve(N_in,fH2O)

t = output['time']
H2 = output['H2']*output['Pres']
CO = output['CO']*output['Pres']
CH4 = output['CH4']*output['Pres']
CO2 = output['CO2']*output['Pres']
H2O = output['H2O']*output['Pres']
NH3 = output['NH3']*output['Pres']
N2 = output['N2']*output['Pres']

Haze=output['dN_Hazedt']
HCN=output['dN_HCNdt']

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
ax2.plot(t/yr/1e6, HCN,'C4:',lw=4)
ax2.plot(t/yr/1e6, Haze,'k:',lw=4)
ax2.set_yscale('log')
ax2.set_ylim(1e7,2e13)

plt.savefig("Fig6_Zahnle2020.pdf",bbox_inches='tight')
plt.show()
