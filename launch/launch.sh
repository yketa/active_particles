#! /bin/bash
#
# This bash shell script submits jobs to a Slurm job scheduler, setting
# environment variables specific to the script.
#
# Jobs are launched as follows:
# [SLURM PARAMETERS] bash launch.sh [COMMAND] [SCRIPT] [ENVIRONMENT VARIABLES]
#
# Slurm parameters
# ----------------
# OUT_DIR : string
#   Error output directory.
#   DEFAULT: active_particles.naming.out_directory
# JOB_NAME : string
#   Job name on Slurm scheduler
#   DEFAULT: script name
# CHAIN : int
#   Begin execution after job with job ID CHAIN has succesfully executed.
#   DEFAULT: None
# PARTITION : string
#   Partition for the resource allocation.
#   DEFAULT: gpu
# GRES : string
#   Generic consumable resources.
#   DEFAULT: gpu:k80:1
# NTASKS : int
#   Maximum ntasks to be invoked on each core.
#   DEFAULT: 1

OUT_DIR=${OUT_DIR-
$(python -c 'from active_particles.naming import out_directory;
print(out_directory)')} # output directory
mkdir -p $OUT_DIR       # create if not existing

COMMAND=$1  # command to execute script
shift
SCRIPT=$1   # script
shift
ENVVAR=$@   # environment variables for script execution

# SUBMIT JOB
sbatch --job-name=${JOB_NAME-${SCRIPT##*/}} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#! /bin/bash
#SBATCH --partition=${PARTITION-gpu}  # partition for the resource allocation
#SBATCH --gres=${GRES-gpu:k80:1}      # generic consumable resources
#SBATCH --error=${OUT_DIR}/%j.out     # standard error output file
#SBATCH --ntasks-per-node=${NTASKS-1} # maximum ntasks to be invoked on each core

export $ENVVAR                                  # setting environment variables
(>&2 printf '%-8s: %s\n' 'COMMAND' '$COMMAND')  # COMMAND in error output file
(>&2 printf '%-8s: %s\n' 'SCRIPT' '$SCRIPT')    # SCRIPT in error output file
(>&2 printf '%-8s: %s\n' 'ENVIRON' '$ENVVAR')   # ENVVAR in error output file
(>&2 echo)

$COMMAND $SCRIPT  # launching script
EOF
