
import re


# note the atoms in CH2CHO doesn't parse properly. So edit zahnle.cti after running
# script to manually put in atoms.

thermo_file = 'thermodata120_2-10-2021.rx'
rx_file = 'titan_205_2-10-2021_modified.rx'

def write_2body(rx,a,b,Ea):
    start = 'reaction( '
    reaction = '"'+rx[0]+' + '+rx[1]+' <=> '+rx[2]+' + '+rx[3]+'",'
    data = '  ['+'%.6e'%(a*(1/300)**b)+', '+'%.3f'%b+', '+'%.2f'%Ea+'])'

    out = ''.join([start,reaction,data])
    
    if any([rrr=="M" for rrr in rx]):
        start = 'three_body_reaction( '
        reaction = '"'+rx[0]+' + '+rx[1]+' <=> '+rx[2]+' + '+rx[3]+'",'
        data = '  ['+'%.6e'%(a*(1/300)**b)+', '+'%.3f'%b+', '+'%.2f'%Ea+'])'

        out = ''.join([start,reaction,data])
    return out
def write_falloff(rx,a,b,Ea0,d,e,Eainf):
    if len(rx)==3:
        start = 'falloff_reaction( '
        reaction = '"'+rx[0]+' + '+rx[1]+' (+ M) <=> '+rx[2]+' (+ M)",\n'
        kf = '         kf  = ['+'%.6e'%(d*(1/300)**e)+', '+'%.3f'%e+', '+'%.2f'%Eainf+'],\n'
        kf0 = '        kf0 = ['+'%.6e'%(a*(1/300)**b)+', '+'%.3f'%b+', '+'%.2f'%Ea0+'])'

        out = ''.join([start,reaction,kf,kf0])
    else:
        start = 'falloff_reaction( '
        reaction = '"'+rx[0]+' + '+rx[1]+' (+ M) <=> '+rx[2]+' + '+rx[3]+' (+ M)",\n'
        kf = '         kf  = ['+'%.6e'%(d*(1/300)**e)+', '+'%.3f'%e+', '+'%.2f'%Eainf+'],\n'
        kf0 = '        kf0 = ['+'%.6e'%(a*(1/300)**b)+', '+'%.3f'%b+', '+'%.2f'%Ea0+'])'

        out = ''.join([start,reaction,kf,kf0])
    return out
    
def return_coefficients(spec):
    fil = open(thermo_file,'r')
    lines = fil.readlines()
    fil.close()
    
    lines1 = []
    for line in lines:
        if line.split()[0] == spec:
            lines1.append(line)
            
    if len(lines)==0:
        raise Exception('species not found')
    
    
    coeff = []
    lb = []
    ub = []
    for line in lines1:
    
        if len(line.split())==4:
            raise Exception('no data')
            pass
        else:
            coeff.append([float(a) for a in line.split()[5:]])
            lb.append(float(line.split()[1])*1000)
            ub.append(float(line.split()[2])*1000)
    
            
    return lb,ub,coeff


def return_coefficients_nodata(spec):
    fil = open(thermo_file,'r')
    lines = fil.readlines()
    fil.close()
    
    for line in lines:
        if line.split()[0]==spec:
            data = line

    dum, xH_j, xS_j, spec_ref = data.split()
    
    for line in lines:
        if line.split()[0]==spec_ref:
            data = line
    
    dum, dum1,dum,xH_i, xS_i = data.split()[:5]
    xH_i = float(xH_i)
    xS_i = float(xS_i)
    xH_j = float(xH_j)
    xS_j = float(xS_j)

    deltaH = (xH_j-xH_i)
    deltaS = (xS_j-xS_i)

    lb,ub,coeff = return_coefficients(spec_ref)
    
    coeff1 = []
    for co in coeff:
        co1 = co.copy()
        co1[5] = co[5] + deltaH
        co1[6] = co[6] + deltaS
        coeff1.append(co1)
    return lb,ub,coeff1,spec_ref
    
    
def find_atoms(a):
    bb = re.findall(r'([A-Z][a-z]?)(\d*)', a)

    bbb = []
    for b in bb:
        bbb.append(list(b))
    bb = bbb

    bbb = []
    jjj = 0
    for i,b1 in enumerate(bb):
        k = 0
        for j,b2 in enumerate(bb):
            if i!=j:
                if b1[0] == b2[0]:
                    k+=1
                    if len(b1[1])==0:
                        b1[1] = '1'
                    if len(b2[1])==0:
                        b2[1]='1'
                    bnew = [b1[0],str(int(b1[1])+int(b2[1]))]

        if k==0:  
            bbb.append(b1)
        elif k==1 and jjj==0:
            bbb.append(bnew)
            jjj+=1
    bb = bbb

    tot = []
    for b in bb:
        if b[0]=='D':
            pass
        else:
            if len(b[1])==0:
                tot.append(b[0]+':1')
            else:
                tot.append(b[0]+':'+b[1])

    tot1 = '  '.join(tot)
    return tot1

def shomate(Tlow,Thigh,coeff):
    return '            Shomate( ['+"{:>8}".format('%.1f'%Tlow)+','+"{:>8}".format('%.1f'%Thigh)+'], [ '+'%.8e'%coeff[0]+', '+'%.8e'%coeff[1]+'\n'+\
           '                    , '+'%.8e'%coeff[2]+', '+'%.8e'%coeff[3]+', '+'%.8e'%coeff[4]+'\n'+\
           '                    , '+'%.8e'%coeff[5]+', '+'%.8e'%coeff[6]+'] )'

def thermo_entry(spec,lb,ub,coeff,note = '"From the NIST database"'):
    
    atoms = find_atoms(spec)
    
    start = 'species(name = "'+spec+'",\n'+\
        '        atoms = " '+atoms+' ",\n'+\
        '        thermo = (\n'
    thermo = []
    for i in range(len(coeff)):
        if i<len(coeff)-1:
            thermo.append(shomate(lb[i],ub[i],coeff[i])+',\n')
        else:
            thermo.append(shomate(lb[i],ub[i],coeff[i])+'\n')

        
    notes = '                 ),\n        note = '+note
    end = ')\n\n'

    entry = [start]+thermo+[notes]+[end]
    return ''.join(entry)
    
    
    
fil = open(rx_file,'r')
lines = fil.readlines()
fil.close()
rxspecies = []
for line in lines:
    for sp in line[:36].split():
        if not any([sp==spp for spp in rxspecies]):
            rxspecies.append(sp)
            
            
fil = open(thermo_file,'r')
lines = fil.readlines()
fil.close()

fil = open('thermo_data.txt','w')

spec = []
for line in lines:
    sp = line.split()[0].strip('AER')
    if not any([sp==spp for spp in spec]):
        if any([sp==spp for spp in rxspecies]):
            spec.append(sp)
            
# print(set(spec).symmetric_difference(set(rxspecies)))
            
index = spec.index('HCS')
for i in range(index):
    lb,ub,coeff = return_coefficients(spec[i])
    # only use two smallest sets of coefficients
    if len(lb)>2:
        lb = lb[0:2]
        ub = ub[0:2]
        coeff = coeff[0:2]
    
    if len(lb)>0:
        res = thermo_entry(spec[i],lb,ub,coeff)
        fil.write(res)
        
# now need to do species without data
for i in range(index,len(spec)):

    lb,ub,coeff,spec_ref = return_coefficients_nodata(spec[i])
    if len(lb)>2:
        lb = lb[0:2]
        ub = ub[0:2]
        coeff = coeff[0:2]
    
    if len(lb)==0:
        sys.exit(spec[i])
    res = thermo_entry(spec[i],lb,ub,coeff,note='"Estimated from thermodynamic data at 298 K and species '+spec_ref+'"')
    fil.write(res)


fil.close() 



# now lets work on reactions
fil = open(rx_file,'r')
lines = fil.readlines()
fil.close()

fil = open('reactions.txt','w')

for i,line in enumerate(lines):
    rx = line[:36].split()
    rates = [float(a) for a in line[36:].split()]
    a,b,c,d,e,f = rates
    if d<=0 and a>0:
        # twobody
        Ea = c
        # now we have
        # a, b, Er
        out = write_2body(rx,a,b,Ea)
        
        fil.write('# Reaction '+str(i+1)+'\n')
        fil.write(out+'\n\n')
        
    else:
        # threebody falloff
        Ea0 = c
        Eainf = f
        out = write_falloff(rx,a,b,Ea0,d,e,Eainf)
        
        fil.write('# Reaction '+str(i+1)+'\n')
        fil.write(out+'\n\n')
        
fil.close()


# combine

header = '''# Generated from titan_205.rx, and thermodata120.rx
# on 4/16/2021. Kevin sent me (Nick Wogan) these data files on 
# 2/10/2021 in an email. This is Kevin Zahnle's reaction mechanism 
# designed for atmospheric chemistry.

'''

units = 'units(length="cm", quantity = "molec", act_energy = "K")\n\n'

temp =[] 
rxspecies.remove('M')

for i in range(int(len(rxspecies)/10)+1):
    if i<7:
        temp.append(' '.join(rxspecies[10*i:10*i+10])+'\n')
    else:
        temp.append(' '.join(rxspecies[10*i:])+'')
try:
    temp.remove('')
except:
    pass
        
    
species = '                       '.join(temp)



ideal_gas = 'ideal_gas(name = "zahnle",\n'+\
            '          elements = "O H C N S",\n'+\
            '          species = """'+species+'""",\n'+\
            '          reactions = "all",\n'+\
            '          initial_state = state(temperature = 300.0, pressure = OneAtm) )\n\n'


temp1 = '''#-------------------------------------------------------------------------------
#  Species data 
#-------------------------------------------------------------------------------
'''

fil = open('thermo_data.txt','r')
lines1 = fil.readlines()
fil.close()


temp2 = '''#-------------------------------------------------------------------------------
#  Reaction data 
#-------------------------------------------------------------------------------
'''

fil = open('reactions.txt','r')
lines2 = fil.readlines()
fil.close()



together = ''.join([header,units,ideal_gas,temp1,''.join(lines1),temp2,''.join(lines2)])

fil = open('zahnle.cti','w')
fil.write(together)
fil.close()

import os
os.remove("thermo_data.txt")
os.remove("reactions.txt")

