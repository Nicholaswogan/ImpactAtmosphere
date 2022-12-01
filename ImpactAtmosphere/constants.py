import cantera as ct

L_H2O = 2307.61404e7 # erg/g
mu_H2O = 18.015 # g/mol
R = 8.31446261815324e7 # erg/(mol*K)
N_avo = ct.avogadro/1e3 # molecules/mol
k_boltz = ct.boltzmann*1e7 # egs/K
yr = 365*24*60*60.0 # s in year
Me = 5.972e24 # mass of Earth in kg

# heat capacity of steam gas
cp_H2O = 1996e4 # erg/(g*K)
T_crit_H2O = 647.096 # K, critical point of H2O

# ~half machine precision tolerance for various things
TOL = 1.0e-9