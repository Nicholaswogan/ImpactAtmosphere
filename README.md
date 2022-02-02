# ImpactAtmosphere
This package contains models of atmospheric evolution after large asteroid impacts on the Hadean Earth.

This repository builds upon [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c) (see the original [source code here](https://zenodo.org/record/3698264#.YCHAqndKhuU)).

# Installation

You need the GNU compiler collection, version >4.9.4 (includes `gfortran`, `gcc`, etc.). If you are using a Mac, I suggest installing it with Homebrew: `brew install gcc`. For other operating systems [follow this GNU installation guide](https://gcc.gnu.org/install/binaries.html).

Next, make a conda environment with the following command

```sh
conda create -n impacts -c cantera -c conda-forge python numpy scipy pyyaml cantera scikit-build cmake ninja matplotlib
```

Activate the environment, navigate to the root directory of this repository and install with

```sh
conda activate impacts
python setup.py install
```

<!-- # Updates
This code builds up on [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c) ([source code here](https://zenodo.org/record/3698264#.YCHAqndKhuU)). Here are the similarities and differences.
- `ImpactAtmosphere.EvolveAtm` and `src/EvolveAtmFort.f90` is a recreation of `photochem_implicit_posted.f` , with the following differences
  - This code has a different hydrogen escape parameterization.
  - This code integrates the ODEs using `scipy.integrate.solve_ivp` (high order/accurate). The version published by [Zahnle et al. (2020)](https://iopscience.iop.org/article/10.3847/PSJ/ab7e2c) used backward Euler (low order method/less accurate).
  - This code self-consistently solves for the UV optical depth at each timestep. The old version of the code used the UV optical depth from the pervious timestep, which sometimes caused the integration to break.

- `ImpactAtmosphere.SteamAtm` is a recreation of `IW_posted.f`, but is fundumentally different

- This code has an additional feature, `EvolveAtm.EvolveAtm.diffuse`, which estimates the HCN mixing ratio as a function of altitude in the atmosphere, given an HCN production rate from photochemistry (in molecules/cm^2/s). -->
