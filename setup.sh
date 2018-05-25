#! /bin/bash
#
# This bash shell script sets up aliases, environment variables and functions
# useful when using active_particles.

export AP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # path to active_particles package

# PYTHON
export PYTHONPATH=$PYTHONPATH:${AP_DIR}/..  # for python to find active_particles

# COMMANDS
alias ap_param="python ${AP_DIR}/param.py"
alias ap_launch="bash ${AP_DIR}/launch/launch.sh"

# SCRIPTS (defined as variables so they can be used with ap_launch)
export AP_CSS="python ${AP_DIR}/analysis/css.py"
export AP_CUU="python ${AP_DIR}/analysis/cuu.py"
export AP_FRAME="python ${AP_DIR}/analysis/frame.py"
export AP_MSD="python ${AP_DIR}/analysis/msd.py"

# FUNCTIONS
. ${AP_DIR}/exponents.sh # translations between floats and litteral expressions
