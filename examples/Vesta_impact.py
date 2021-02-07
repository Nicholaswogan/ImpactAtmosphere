import numpy as np
import cantera as ct
from ImpactAtmosphere import SteamAtm
from matplotlib import pyplot as plt

gas = ct.Solution('gri30.cti') # use Gri-Mech 3.0
N_H2O_ocean = 15.0e3 # moles/cm2 (15e3 is 1 ocean)
N_CO2 = 23.  # moles/cm2 (23 is "one bar" of CO2)
N_N2  = 36. # moles/cm2 (36 is "one bar" of N2)
M_i = 2.589e23 # impactor mass in grams (2.589e23 g = 500 km, about Vesta)
stm = SteamAtm(gas)
solution = stm.impact(N_H2O_ocean,N_CO2,N_N2,M_i)

# plot results
plt.rcParams.update({'font.size': 15})
fig,ax = plt.subplots(1,1,figsize=[9,5])
yr = 365*24*60*60
species = ['CH4','CO','CO2','NH3']
for i,spec in enumerate(species):
    ax.plot(solution['time']/yr,solution[spec]*solution['Ntot'],'C'+str(i)+'-',label=spec)
ax.set_yscale('log')
ax.legend(ncol=2)
ax1 = ax.twinx()
ax1.plot(solution['time']/yr,solution['Tsurf'],'r-',lw=2,label='Temperature')
ax1.legend()
ax.set_ylabel('Column abundance (mol/cm$^2$)')
ax.set_xlabel('Time (years)')
ax1.set_ylabel('Temperature (K)')
plt.show()
