# EvolveAtm
This program is a 0-D photochemical model designed to simulate the Hadean Earth atmosphere after an asteroid impact.

The code was originally developed by Kevin Zahnle and was described in [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c). This repository is an updated version of the code.

# Installation and Use
The code is written in Fortran, but this repository also contains a Python wrapper to the Fortran code. Following are instructions for using the Fortran or Python versions.

## Python
To install python version, download this repository and navigate a terminal to the directory `EvolveAtm`, then install the python package with pip:
```bash
pip install .
```

To check that your installation worked propoerly, run the example in `EvolveAtm/Example`.

## Fortran
The file `photochem_implicit_nickv2.2.f`, is the main fortran program. Using this program all by itself requires generating some input files (see the fortran code for details).

The fortran code can also be run using the bash script `impact_atm_photochem.sh`. This script first runs the program `IW.f`. This program models the cooling of a steamy atmosphere after an impact (See section 3 in [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c) for model details). The script then uses output from `IW.f` as input for `photochem_implicit_nickv2.2.f`. Many outputfiles are made.
