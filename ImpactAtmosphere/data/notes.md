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
CH3CHO  H       CH3CO   H2              6.6E-11   0.0    2100.    0.0       0.0   0.0 # duplicate
C2H4    CH3     C3H6    H               1.2e-11   0.0    5600.    0.0       0.0   0.0 # duplicate. Also updated rate.
N       C2H3    NH      C2H2            1.20E-11  0.0     0.0     0.0       0.0   0.0 # duplicate
H       HNCO    H2      NCO             1.9E-12   1.66   7000.    0.0       0.0   0.0 # duplicate
H2CN    N       NH      C2H2            1.2e-11   0.0      0.0    0.0       0.0   0.0 # bad mass balance
```
