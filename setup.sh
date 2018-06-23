#! /bin/bash
#
# This bash shell script adds the directory containing active_particles to the
# Python path, and sets up aliases, environment variables and functions useful
# when using active_particles.

export AP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # path to active_particles package
alias ap_python=$(. activate active_particles_env; which python)  # python executable to use with active_particles
alias ap_update="( cd $AP_DIR ; git pull )"                       # alias to update active_particles git repository

# PYTHON
export PYTHONPATH=$PYTHONPATH:${AP_DIR}/..  # for python to find active_particles

# COMMANDS
alias ap_param="ap_python ${AP_DIR}/param.py"
alias ap_launch="bash ${AP_DIR}/launch/launch.sh"

# SCRIPTS (defined as variables so they can be used with ap_launch)
export AP_CSS="ap_python ${AP_DIR}/analysis/css.py"
export AP_CTT="ap_python ${AP_DIR}/analysis/ctt.py"
export AP_CUU="ap_python ${AP_DIR}/analysis/cuu.py"
export AP_FRAME="ap_python ${AP_DIR}/analysis/frame.py"
export AP_MSD="ap_python ${AP_DIR}/analysis/msd.py"
export AP_VARN="ap_python ${AP_DIR}/analysis/varn.py"

# PLOTS
alias ap_c44="ap_python ${AP_DIR}/plot/c44.py"
alias ap_pphiloc="ap_python ${AP_DIR}/plot/pphiloc.py"

# FUNCTIONS
. ${AP_DIR}/exponents.sh # translations between floats and litteral expressions
