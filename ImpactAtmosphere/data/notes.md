# Kevin Zahnle's reaction mechanism

## 4/16/2021
Generated new .cti and .yaml file formats of Kevin Zahnle's reaction mechanism. I used reaction tables and thermodynamic data that he sent me via email on 2/10/2021 (`titan_205_2-10-2021.rx`, and `thermodata120_2-10-2021.rx`). However, before converting i deleted the following reactions. The file that does not contain these reactions is `titan_205_2-10-2021_modified.rx`.
```
H       H2COH   CH3     OH              5.0e-11   0.0     0.0     0.0       0.0   0.0 # duplicate
H       HNCO    H2      NCO             1.9E-12   1.66   7000.    0.0       0.0   0.0 # duplicate
NO2     OH      HO2     NO              3.0E-11   0.0    3360.    0.0       0.0   0.0 # duplicate
C2H3    N       NH      C2H2            1.2e-11   0.0      0.0    0.0       0.0   0.0 # duplicate  
H2CN    N       NH      C2H2            1.2e-11   0.0      0.0    0.0       0.0   0.0 # bad mass balance
```

Additionally, Cantera did not allow 3 temperature ranges of `Shomate` thermodynamic data. Therefore, I did not use Kevin's highest temperature data. I'm guessing he got this data from the NASA polynomials, then converted it to Shomate format.

Converting to .cti with `python zahnle2cti.py`. However this script can't parse the atoms in CH2CHO, so after conversion, manually correct the number of atoms in this molecule. Convert to .yaml with `python -m cantera.cti2yaml zahnle.cti`.

To run the validation of these reaction mechanisms run `python validation.py`.

## 6/6/2021
Kevin sent me `earth_125_5-26-2021.rx` via email on 5/26/2021. I deleted all the photolysis reactions, and several other reactions in this file. The file that does not contain these reactions is `earth_125_5-26-2021_modified.rx`. Here are the non-photolysis reactions I deleted, and the reasons
```
H2      CH3CO   CH3CHO  H               2.18E-13  1.82   8860.    0.0       0.0   0.0 # duplicate
C3H6    H       C2H4    CH3             2.2e-12   1.5    1010.    0.0       0.0   0.0 # duplicate
N       C2H3    NH      C2H2            1.20E-11  0.0     0.0     0.0       0.0   0.0 # duplicate
H       HNCO    H2      NCO             1.9E-12   1.66   7000.    0.0       0.0   0.0 # duplicate
H2CN    N       NH      C2H2            1.2e-11   0.0      0.0    0.0       0.0   0.0 # bad mass balance
```

I also updated a few reaction rates
```
CH3CHO  H       CH3CO   H2              6.6E-11   0.0    2100.    0.0       0.0   0.0 # Old rate
CH3CHO  H       CH3CO   H2              5.27e-13  2.58    615.    0.0       0.0   0.0 # New rate

C2H4    CH3     C3H6    H               1.2e-11   0.0    5600.    0.0       0.0   0.0 # Old rate
C2H4    CH3     C3H6    H               6.9e-12   0.0    5600.    0.0       0.0   0.0 # New rate
```

# 9/7/21

I made two versions of `zahnle_earth.yaml`. One version has S8 and S8L and the other does not. S8 and S8L are very fast species which cause integration problems. Default is no S8.
