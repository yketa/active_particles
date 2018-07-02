#! /bin/bash
#
# This bash shell script adds the directory containing active_particles to the
# Python path, and sets up aliases, environment variables and functions useful
# when using active_particles.

export AP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # path to active_particles package
export AP_PYTHON=$(. activate active_particles_env; which python) # python executable to use with active_particles
alias ap_python="$AP_PYTHON"
alias ap_update="( cd $AP_DIR ; git pull ; . setup.sh )"          # alias to update active_particles git repository

# PYTHON
export PYTHONPATH=$PYTHONPATH:${AP_DIR}/..  # for python to find active_particles

# COMMANDS
alias ap_param="$AP_PYTHON ${AP_DIR}/param.py"
alias ap_launch="bash ${AP_DIR}/launch/launch.sh"

# SCRIPTS (defined as variables so they can be used with ap_launch)
export AP_CSS="$AP_PYTHON ${AP_DIR}/analysis/css.py"
export AP_CTT="$AP_PYTHON ${AP_DIR}/analysis/ctt.py"
export AP_CUU="$AP_PYTHON ${AP_DIR}/analysis/cuu.py"
export AP_FRAME="$AP_PYTHON ${AP_DIR}/analysis/frame.py"
export AP_MSD="$AP_PYTHON ${AP_DIR}/analysis/msd.py"
export AP_VARN="$AP_PYTHON ${AP_DIR}/analysis/varn.py"

# PLOTS
alias ap_c44="$AP_PYTHON ${AP_DIR}/plot/c44.py"
alias ap_chi_msd="$AP_PYTHON ${AP_DIR}/plot/chi_msd.py"
alias ap_pphiloc="$AP_PYTHON ${AP_DIR}/plot/pphiloc.py"

# FUNCTIONS
. ${AP_DIR}/exponents.sh # translations between floats and litteral expressions
