![Active Brownian disks](https://github.com/yketa/UBC_2018_Wiki/raw/master/Images/header.png)

# Simple model of active particles
### Yann-Edwin Keta, University of British Columbia, 2018

This repository contains scripts, for simulation and analysis purposes, developed for a research project, detailed in __[this wiki](https://yketa.github.io/UBC_2018_Wiki)__. It is organised as a Python importable package.

Simulation scripts rely on the [HOOMD-blue](https://glotzerlab.engin.umich.edu/hoomd-blue/) simulation toolkit to perform molecular dynamics simulation of active particles.

## Requirements

* __active_particles__ was tested on Mac OS X and Linux, thus the installation script should work on these OS only.
* An environment is set specifically for __active_particles__ with [conda](https://conda.io/miniconda.html), which thus has to be installed.

## Installation

Clone this repository and source the [`install.sh`](https://github.com/yketa/active_particles/blob/master/install.sh) file.
```bash
git clone https://github.com/yketa/active_particles
source active_particles/install.sh
```

[`install.sh`](https://github.com/yketa/active_particles/blob/master/install.sh) creates a specific conda environment, sources [`setup.sh`](https://github.com/yketa/active_particles/blob/master/setup.sh) — which adds the directory containing __active_particles__ to the Python path, and sets up aliases, environment variables and functions to be used in the command line — and writes to `~/.bash_profile` so the latter is done at each login.

__Remark:__ This is not an installation _stricto sensu_, all scripts can be used as they come without running [`install.sh`](https://github.com/yketa/active_particles/blob/master/install.sh). However, running it then using the defined commands and environment variables makes sure everything works as expected.

## Useful aliases

* `ap_python`: Python executable installed in the __active_particles__ conda environment.
* `ap_mprof`: `mprof` command of [`memory_profiler`](https://github.com/pythonprofilers/memory_profiler), useful for monitoring memory consumption of processes, installed in the __active_particles__ conda environment.
* `ap_update`: Pull git repository and source [`setup.sh`](https://github.com/yketa/active_particles/blob/master/setup.sh).

## Notable scripts

* [`setup.sh`](https://github.com/yketa/active_particles/blob/master/setup.sh) sets up aliases, environment variables and functions to be used in the command line.
* [`naming.py`](https://github.com/yketa/active_particles/blob/master/naming.py) contains all the standards for naming all data files and directories.

## Usage

All modules (or almost) contain a detailed documentation in their header, accessible through Python with the built-in function `help(module_name)`. Scripts' headers should contain a detailed list of parameters to set as environment variables when executing them.

For example, to launch the script [`analysis/cuu.py`](https://github.com/yketa/active_particles/blob/master/analysis/cuu.py) — which computes displacement correlations — one can use the environment variable `AP_CUU` set by [`setup.sh`](https://github.com/yketa/active_particles/blob/master/setup.sh). With parameters:

* compute displacement correlations = True,
* lag time = 10,
* numbers of intervals of time considered = 10,

and all other parameters set to default, the program is launched the following way
```bash
yketa@yketa: ~/active_particles_data/Dk8000_Vj1000_Rg2000_Nq1000_Ll0000 $ COMPUTE=True TIME=10 INTERVAL_MAXIMUM=10 $AP_CUU
Execution time: 0:01:26.969770
```
and will output its execution time. Launching the program to plot the calculated data is done similarly
```bash
yketa@yketa: ~/active_particles_data/Dk8000_Vj1000_Rg2000_Nq1000_Ll0000 $ SHOW=True PLOT=True R_MAX=50 AXIS=LINLOG TIME=10 INTERVAL_MAXIMUM=10 $AP_CUU
```
and will make plotting windows appear.

![Example of Cuu plot](https://github.com/yketa/active_particles/raw/master/docs/screenshot_cuu.png)

