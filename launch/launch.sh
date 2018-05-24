#! /bin/bash
#
# This bash shell script submits jobs to a Slurm job scheduler, setting
# environment variables specific to the script.
#
# Jobs are launched as follows:
# [SLURM PARAMETERS] bash launch.sh [COMMAND] [SCRIPT] [ENVIRONMENT VARIABLES]

OUT_DIR=$(python -c 'from active_particles.naming import out_directory;
print(out_directory)')  # output directory
mkdir -p $OUT_DIR       # create if not existing

COMMAND=$1  # command to execute script
shift
SCRIPT=$1   # script
shift
ENVVAR=$@   # environment variables for script execution

# SUBMIT JOB
sbatch --job-name=${SCRIPT##*/} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#! /bin/bash
#SBATCH --partition=${PARTITION-gpu}  # partition for the resource allocation
#SBATCH --gres=${GRES-gpu:k80:1}      # generic consumable resources
#SBATCH --error=${OUT_DIR}/%j.out     # standard error output file
#SBATCH --ntasks-per-node=${NTASKS-1} # maximum ntasks to be invoked on each core

export $ENVVAR                            # setting environment variables
(>&2 echo COMMAND: $COMMAND)              # COMMAND in error output file
(>&2 echo SCRIPT: $SCRIPT)                # SCRIPT in error output file
(>&2 echo ENVIRONMENT VARIABLES: $ENVVAR) # ENVVAR in error output file

$COMMAND $SCRIPT  # launching script
EOF
