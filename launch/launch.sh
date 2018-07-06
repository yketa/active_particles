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
# CHAIN : int
#   Begin execution after job with job ID CHAIN has succesfully executed.
#   DEFAULT: None
# JOB_NAME : string
#   Job name on Slurm scheduler
#   DEFAULT: script name
# PARTITION : string
#   Partition for the resource allocation.
#   DEFAULT: vis
# GRES : string
#   Generic consumable resources.
#   DEFAULT:
# OUT_FILE : string
#   Standard output file.
#   DEFAULT: /dev/null (output file is supposed to be managed in script itself)
# NTASKS : int
#   Maximum ntasks to be invoked on each core.
#   DEFAULT: 1

OUT_DIR=${OUT_DIR-$($AP_PYTHON -c 'from active_particles.naming import out_directory; print(out_directory)')} # output directory
mkdir -p $OUT_DIR                                                                                             # create if not existing

COMMAND=$1  # command to execute script
if [[ -z "$COMMAND" ]]; then
  echo 'No command submitted.'
  exit 0
fi
shift
SCRIPT=$1   # script
if [[ -z "$SCRIPT" ]]; then
  echo 'No script submitted.'
  exit 0
fi
shift
ENVVAR=$@   # environment variables for script execution

if [[ ! -z "$DATA" ]]; then # data directory name submitted
  SIM_DIR=$($AP_PYTHON -c 'from active_particles.naming import sim_directory; print(sim_directory)')
  ENVVAR="DATA_DIRECTORY=${SIM_DIR}/${DATA} $ENVVAR"
fi

# SUBMIT JOB
sbatch ${CHAIN:+-d afterok:$CHAIN} <<EOF
#! /bin/bash
#SBATCH --job-name=${JOB_NAME-${SCRIPT##*/}}  # job name
#SBATCH --partition=${PARTITION-vis}          # partition for the resource allocation
#SBATCH --gres=${GRES-}                       # generic consumable resources
#SBATCH --output=${OUT_FILE-/dev/null}        # standard output file
#SBATCH --error=${OUT_DIR}/%j.out             # standard error output file
#SBATCH --ntasks-per-node=${NTASKS-1}         # maximum ntasks to be invoked on each core

export $ENVVAR  # setting environment variables

(>&2 printf '%-17s: %s\n' 'SUBMIT DIRECTORY' '$(pwd)')  # submit directory in error output file
(>&2 printf '%-17s: %s\n' 'COMMAND' '$COMMAND')         # COMMAND in error output file
(>&2 printf '%-17s: %s\n' 'SCRIPT' '$SCRIPT')           # SCRIPT in error output file
(>&2 printf '%-17s: %s\n' 'ENVIRON' '$ENVVAR')          # ENVVAR in error output file
(>&2 echo)

$COMMAND $SCRIPT  # launching script
EOF
