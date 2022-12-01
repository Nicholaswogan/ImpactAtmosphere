from scipy import interpolate
import numpy as np
import os
from . import constants as const

data_dir = os.path.dirname(os.path.realpath(__file__))+'/data/'

def read_table_C3():
    col_labels = ['Mproj_Mearth','vimp_vesc','theta','Eimp','M_Fe_proj','X_Fe_interior','X_Fe_surf','X_Fe_atmos','X_Fe_disk','X_Fe_ejec']
    out = {}
    for col in col_labels:
        out[col] = []

    with open(data_dir+'Citron_2022_C3.txt','r') as f:
        lines = f.readlines()    
        for line in lines[5:53]:
            tmp = [float(a.replace(' ','').replace('x10^','e')) for a in line.strip().split('\t')]

            for i,col in enumerate(col_labels):
                out[col].append(tmp[i])

    v = out['vimp_vesc']
    v = list(set(v))
    theta = out['theta']
    theta = list(set(theta))
    M = list(set(out['Mproj_Mearth']))

    out1 = {}

    for vv in v:
        out1[vv] = {}
        for thetaa in theta:
            out1[vv][thetaa] = {}
            for MM in M:

                tmp = {}
                for i in range(len(out['theta'])):
                    if out['theta'][i] == thetaa and out['vimp_vesc'][i] == vv and out['Mproj_Mearth'][i] == MM:
                        for col in col_labels[3:]:
                            tmp[col] = out[col][i]
                        break
                out1[vv][thetaa][MM] = tmp
            
    return out, out1

def read_table_C2():
    col_labels = ['Mproj_Mearth','vimp_vesc','theta','Eimp','M_melt','M_scf','M_scf_atmos','M_vapor','M_atmos','M_disk']
    out = {}
    for col in col_labels:
        out[col] = []

    with open(data_dir+'Citron_2022_C2.txt','r') as f:
        lines = f.readlines()    
        for line in lines[5:53]:
            tmp = [float(a.replace(' ','').replace('x10^','e')) for a in line.strip().split('\t')]

            for i,col in enumerate(col_labels):
                out[col].append(tmp[i])

    v = out['vimp_vesc']
    v = list(set(v))
    theta = out['theta']
    theta = list(set(theta))
    M = list(set(out['Mproj_Mearth']))

    out1 = {}

    for vv in v:
        out1[vv] = {}
        for thetaa in theta:
            out1[vv][thetaa] = {}
            for MM in M:

                tmp = {}
                for i in range(len(out['theta'])):
                    if out['theta'][i] == thetaa and out['vimp_vesc'][i] == vv and out['Mproj_Mearth'][i] == MM:
                        for col in col_labels[3:]:
                            tmp[col] = out[col][i]
                        break
                out1[vv][thetaa][MM] = tmp
            
    return out, out1

def read_table_C1():
    col_labels = ['Mproj_Mearth','vimp_vesc','theta','Eimp','X_E_interior','X_E_surf','X_E_atmos','X_E_disk','X_E_esc','X_IE_surf','X_IE_atmos']
    out = {}
    for col in col_labels:
        out[col] = []

    with open('ImpactAtmosphere/data/'+'Citron_2022_C1.txt','r') as f:
        lines = f.readlines()    
        for line in lines[5:53]:
            tmp = [float(a.replace(' ','').replace('x10^','e')) for a in line.strip().split('\t')]

            for i,col in enumerate(col_labels):
                out[col].append(tmp[i])

    v = out['vimp_vesc']
    v = list(set(v))
    theta = out['theta']
    theta = list(set(theta))
    M = list(set(out['Mproj_Mearth']))

    out1 = {}

    for vv in v:
        out1[vv] = {}
        for thetaa in theta:
            out1[vv][thetaa] = {}
            for MM in M:

                tmp = {}
                for i in range(len(out['theta'])):
                    if out['theta'][i] == thetaa and out['vimp_vesc'][i] == vv and out['Mproj_Mearth'][i] == MM:
                        for col in col_labels[3:]:
                            tmp[col] = out[col][i]
                        break
                out1[vv][thetaa][MM] = tmp
            
    return out, out1

# nominal interpolators. We assume 45 degree impact angle, and v_esc = 2.0.
# We linearly interpolate iron and energy results. We linearly interpolate
# rock results, except mass is log10 space.

def nominal_iron_interpolator():
    v = 2 # 
    theta = 45
    col = 'X_Fe_atmos'

    _, iron1 = read_table_C3()
    tmp = iron1[v][theta]
    Mi = []
    XX = []
    for key in tmp:
        Mi.append(key*const.Me*1.0e3) # grams
        XX.append(tmp[key][col])

    iron = interpolate.interp1d(Mi, XX, fill_value='extrapolate')
    return iron

def nominal_rock_interpolator():
    v = 2 # 
    theta = 45

    _, rock1 = read_table_C2()
    tmp = rock1[v][theta]
    Mi = []
    XX = []
    for key in tmp:
        Mi.append(key*const.Me*1.0e3) # grams
        XX.append((tmp[key]['M_melt'] + tmp[key]['M_scf'] + tmp[key]['M_vapor'])*1e3)

    rock = interpolate.interp1d(np.log10(Mi), np.log10(XX), fill_value='extrapolate')
    return rock

def nominal_energy_interpolator():
    v = 2 # 
    theta = 45
    col = 'X_IE_atmos'

    _, energy1 = read_table_C1()
    tmp = energy1[v][theta]
    Mi = []
    XX = []
    for key in tmp:
        Mi.append(key*const.Me*1.0e3) # grams
        XX.append(tmp[key][col])

    rock = interpolate.interp1d(Mi, XX, fill_value='extrapolate')
    return rock


