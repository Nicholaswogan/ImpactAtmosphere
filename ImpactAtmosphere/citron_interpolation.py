from scipy import interpolate
import os

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

def make_interpolator(out1, v, theta, k):
    Me = 5.972e24 # kg in Earth
    tmp = out1[v][theta]
    Mi = []
    XX = []
    for key in tmp:
        Mi.append(key*Me*1.0e3) # grams
        XX.append(tmp[key][k])
    
    f = interpolate.interp1d(Mi, XX, fill_value='extrapolate')
    return f

def make_iron_interpolator(v, theta, col):
    _, iron = read_table_C3()
    Me = 5.972e24 # kg in Earth
    tmp = iron[v][theta]
    Mi = []
    XX = []
    for key in tmp:
        Mi.append(key*Me*1.0e3) # grams
        XX.append(tmp[key][col])
    
    f = interpolate.interp1d(Mi, XX, fill_value='extrapolate')
    return f

def make_rock_interpolator(v, theta, col):
    _, rock = read_table_C2()
    Me = 5.972e24 # kg in Earth
    tmp = rock[v][theta]
    Mi = []
    XX = []
    for key in tmp:
        Mi.append(key*Me*1.0e3) # grams
        XX.append(tmp[key][col])
    
    f = interpolate.interp1d(Mi, XX, fill_value='extrapolate')
    return f



