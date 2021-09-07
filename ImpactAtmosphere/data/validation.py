import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

validation1 = '''Here, I find chemical equilibrium of a gas with
equal parts C, O, H, N with both the zahnle.cti and gri30.cti reaction
mechanisms at T = 1000 K and P = 1 bar. The table shows the difference between
reaction mechanisms are given equal parts C, O, H, N, then is
equilibrated (chemical equilibrium). This shows the differences 
between the two equilibriums
'''
print(validation1)
gas1 = ct.Solution('zahnle_earth.yaml')
gas2 = ct.Solution('gri30.cti')
gas1.TPX = 1000,100000, {'C':1,'O':1,'H':1,'N':1}
gas2.TPX = 1000,100000, {'C':1,'O':1,'H':1,'N':1}
gas1.equilibrate('TP')
gas2.equilibrate('TP')

print("{:>8}".format('species'),"{:>8}".format('ratio'),"{:>8}".format('zahnle'),"{:>8}".format('gri30'))
for key1 in gas2.mole_fraction_dict().keys():
    for key2 in gas1.mole_fraction_dict().keys():
        if key1==key2:
            val1 = gas1.mole_fraction_dict()[key1]
            val2 = gas2.mole_fraction_dict()[key1]
            print("{:>8}".format(key1),'%.2e'%(val1/val2),'%.2e'%val1,'%.2e'%val2)
            
validation2 = '''Here, I find chemical equilibrium of a gas with
equal parts C, O, H, N with both the zahnle.cti and gri30.cti reaction
mechanisms at T = 300 K and P = 1 bar. I then compute the forward
and reverse reaction rates of two body reactions to make sure none of them
are higher than the theoretical maximum limit. See the plot.
'''
print(validation2)
gas = ct.Solution('zahnle_earth.yaml')
P = 100000
T = 300
gas.TPX = T, P, 'CH4:1.0, O2:0.1, N:1, S:1'
gas.equilibrate('TP')
# gas()
# find only 2body reactions
index = []
for i in range(len(gas.forward_rate_constants)):
    if gas.reaction_type(i)==1:
        index.append(i)
        

gasg = ct.Solution('gri30.cti')
gasg.TPX = T, P, 'CH4:1.0, O2:0.1, N:1'
gasg.equilibrate('TP')
# gas()
indexg = []
for i in range(len(gasg.forward_rate_constants)):
    if gasg.reaction_type(i)==1:
        indexg.append(i)
    
plt.rcParams.update({'font.size': 15})
f, [ax,ax1] = plt.subplots(1,2,figsize=[14,5])

ax1.semilogy(index,gas.forward_rate_constants[index]/(ct.avogadro/1e6), 'C0.', label='forward',alpha=1)
ax1.semilogy(index,gas.reverse_rate_constants[index]/(ct.avogadro/1e6), 'C1.', label='reverse',alpha=1)
ax1.set_ylim(1e-30,15.09364569147317)

ax.semilogy(indexg,gasg.forward_rate_constants[indexg]/(ct.avogadro/1e6), 'C0.', label='forward')
ax.semilogy(indexg,gasg.reverse_rate_constants[indexg]/(ct.avogadro/1e6), 'C1.', label='reverse')
ax.set_ylim(1e-30,15.09364569147317)

limit = 2.9e12/(ct.avogadro/1e6)
ax.axhline(limit,ls='--',color='k')
ax.text(0,limit+3*limit,'Theoretical upper limit')
ax1.axhline(limit,ls='--',color='k')
# ax1.text(0,limit+3*limit,'Theoretical upper limit')

ax.set_ylabel('Rate constant\n(cm$^3$ s molecules$^{-1}$)')
ax1.set_ylabel('Rate constant\n(cm$^3$ s molecules$^{-1}$)')
ax.set_xlabel('2-body reaction number')
ax1.set_xlabel('2-body reaction number')
ax.set_title('Gri-Mech 3.0\nT = '+'%.1f'%T+' K')
ax1.set_title('zahnle.cti\nT = '+'%.1f'%T+' K')

# ax.axis(ymin=1e-30)
ax.legend()
plt.subplots_adjust(wspace=.3)
# plt.savefig("Gri-Mech3.0_vs_titan_205.rx_updated.pdf",bbox_inches='tight')
plt.show()
            
            
