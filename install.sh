#! /bin/bash

AP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"  # path to active_particles package

printf "
active_particles
----------------
Yann-Edwin Keta, University of British Columbia, 2018

****************

This program will install the dependencies of the active_particles package and set it up.

"

# CREATE ACTIVE_PARTICLES CONDA ENVIRONMENT

printf "Creating conda environment: active_particles_env\n"
conda env create --force -f ${AP_DIR}/environment.yml

# SET UP ACTIVE_PARTICLES

printf "\nSetting up package.\n"

# add setup command to bash profile
printf "
# ACTIVE_PARTICLES (added with active_particles installer on $( date +%Y-%m-%d' at '%H:%M:%S ))
. ${AP_DIR}/setup.sh\n\n" >> ~/.bash_profile
# set up for this session
. ${AP_DIR}/setup.sh

# DONE

echo "Installation complete."

