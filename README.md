![Active Brownian disks](https://github.com/yketa/UBC_2018_Wiki/raw/master/Images/header.png)

# Simple model of active particles
### Yann-Edwin Keta, University of British Columbia, 2018

This repository contains scripts, for simulation and analysis purposes, developed for a research project, detailed in __[this wiki](https://yketa.github.io/UBC_2018_Wiki)__. It is organised as a Python importable package.

Simulation scripts rely on the [HOOMD-blue](https://glotzerlab.engin.umich.edu/hoomd-blue/) simulation toolkit to perform molecular dynamics simulation of active particles.

All scripts (or almost) contain a detailed documentation in their header, accessible through Python with the built-in function `help(module_name)`.

## Requirements

* __active_particles__ was tested on Mac OS X and Linux, thus the installation script should work on these OS only.
* An environment is set specifically for __active_particles__ with [conda](https://conda.io/miniconda.html), which thus have to be installed.

## Installation

Clone this repository and source the [`install.sh`](https://github.com/yketa/active_particles/blob/master/install.sh) file.
```bash
git clone https://github.com/yketa/active_particles
source active_particles/install.sh
```

[`install.sh`](https://github.com/yketa/active_particles/blob/master/install.sh) creates a specific conda environment, sources [`setup.sh`](https://github.com/yketa/active_particles/blob/master/setup.sh) — which adds the directory containing __active_particles__ to the Python path, and sets up commands, environment variables and functions — and write to `~/.bash_profile` so the latter is done at each login.

__Remark:__ This is not an installation _stricto sensu_, all scripts can be used as they come without running [`install.sh`](https://github.com/yketa/active_particles/blob/master/install.sh). However, running it then using the defined commands and environment variables makes sure everything works as expected.

## Useful commands

* `ap_python`: Python executable in the __active_particles__ conda environment.
* `ap_update`: Pull git repository.

## Important modules

* [`active_particles.naming`](https://github.com/yketa/active_particles/blob/master/naming.py) contains all the standards for naming all data files and directories.
