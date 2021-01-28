# EvolveAtm
This program is a 0-D photochemical model designed to simulate the Hadean Earth atmosphere after an asteroid impact.

The code was originally developed by Kevin Zahnle and was described in [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c). This repository is an updated version of the code.

# Installation and Use
To install this python version download this repository and navigate a terminal to the directory `EvolveAtm`, then install the python package with pip:

```bash
python -m pip install .
```

# Updates
The code is not exactly the same as what was used for [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c). There are a few differences:
- This code has a different hydrogen escape parameterization.
- This code integrates the ODEs using `scipy.integrate.solve_ivp` (high order/accurate). The version published by [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c) used backward Euler (low order method/less accurate).
- This code self-consistently solves for the UV optical depth at each timestep. The old version of the code used the UV optical depth from the pervious timestep, which sometimes caused the integration to break.
- This code has an additional feature, `EvolveAtm.EvolveAtm.diffuse`, which estimates the HCN mixing ratio as a function of altitude in the atmosphere, given an HCN production rate from photochemistry (in molecules/cm^2/s).
