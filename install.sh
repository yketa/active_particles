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

echo "Creating conda environment: active_particles_env"
command -v conda >/dev/null 2>&1 || { echo >&2 "conda not installed. Please visit https://conda.io/miniconda.html. Aborting."; exit 0; }
conda env create --force -f ${AP_DIR}/environment.yml

# SET UP ACTIVE_PARTICLES

echo "Setting up package."

# add setup command to bash profile
printf "
# ACTIVE_PARTICLES (added with active_particles installer on $( date +%Y-%m-%d' at '%H:%M:%S ))
. ${AP_DIR}/setup.sh\n\n" >> ~/.bash_profile
# set up for this session
. ${AP_DIR}/setup.sh

# CREATING SIMULATION DIRECTORY

echo "Creating simulation directory in home directory."
mkdir -p $(ap_python -c "from active_particles.naming import sim_directory; print(sim_directory)")

# DONE

echo "Installation complete."

