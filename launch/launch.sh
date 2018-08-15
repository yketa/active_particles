#! /bin/bash
#
# Submit job to a Slurm job scheduler.

usage(){  # help menu

underline=$(tput smul)
bold=$(tput bold)
normal=$(tput sgr0)

less <<< "Submit job to a Slurm job scheduler.

${bold}SYNOPSIS${normal}

  launch.sh [options] [script]

${bold}OPTIONS${normal}

  ${underline}without argument${normal}
  -h    Display help menu.
  -r    Run with mprof memory profiler.
        (see https://github.com/pythonprofilers/memory_profiler)

  ${underline}with argument${normal}
  -j    Job name on Slurm scheduler.
        DEFAULT: script name after last '/'
  -o    Error output directory.
        DEFAULT: active_particles.naming.out_directory
  -f    Standard output file.
        DEFAULT: /dev/null (output file is supposed to be managed in script
        itself)
  -c    Begin execution after job with this ID has succesfully executed.
        DEFAULT: (not specified)
  -p    Partition for the resource allocation.
        DEFAULT: vis
  -g    Generic consumable resources.
        DEFAULT:
  -n    Maximum ntasks to be invoked on each core.
        DEFAULT: 1
  -m    Real memory required per node.
        DEFAULT: (not specified)
  -d    Data directory.
        DEFAULT: (not specified)
"

}

while getopts "hj:o:f:c:p:g:n:m:d:r" OPTION; do
  case $OPTION in

    h)  # help menu
      usage ; exit ;;
    j)  # job name
      JOB_NAME=$OPTARG ;;
    o)  # error output directory
      OUT_DIR=$OPTARG ;;
    f)  # error output file
      OUT_FILE=$OPTARG ;;
    c)  # chained job
      CHAIN=$OPTARG ;;
    p)  # partition
      PARTITION=$OPTARG ;;
    g)  # generic consumable resources
      GRES=$OPTARG ;;
    n)  # maximum ntasks
      NTASKS=$OPTARG ;;
    m)  # real memory
      MEM=$OPTARG ;;
    d)  # data directory
      SIM_DIR=$($AP_PYTHON -c 'from active_particles.naming import sim_directory; print(sim_directory)')
      DATA="${SIM_DIR}/$OPTARG"
      ;;
    r)  # run with mprof
      MPROF= ;;

  esac
done
shift $(expr $OPTIND - 1)

if [[ -z "$@" ]]; then
  echo 'No script submitted.'
  exit 0
fi

SCRIPT=( ${DATA:+DATA_DIRECTORY=$DATA } $@ )  # script to execute
if [[ ! -z "${MPROF+MPROF}" ]]; then          # run with mprof
  for ((i=0; i<${#SCRIPT[@]}; i++)); do
    if [[ ! "${SCRIPT[$i]}" =~ "=" ]]; then
      break
    fi
  done
  SCRIPT=( ${SCRIPT[@]::$i} $AP_MPROF run ${SCRIPT[@]:$i} )
fi
SCRIPT=${SCRIPT[@]}

# JOB PARAMETERS

OUT_DIR=${OUT_DIR-$($AP_PYTHON -c 'from active_particles.naming import out_directory; print(out_directory)')} # output directory
mkdir -p $OUT_DIR                                                                                             # create if not existing

JOB_NAME=${JOB_NAME-${SCRIPT##*/}}  # job name
PARTITION=${PARTITION-vis}          # partition for the resource allocation
GRES=${GRES-}                       # generic consumable resources
OUT_FILE=${OUT_FILE-/dev/null}      # standard output file
NTASKS=${NTASKS-1}                  # maximum ntasks to be invoked on each core
MEM=${MEM-}                         # real memory required per node

# SUBMIT JOB

sbatch ${CHAIN:+-d afterok:$CHAIN} <<EOF
#! /bin/bash
#SBATCH --job-name='$JOB_NAME'
#SBATCH --partition=$PARTITION
#SBATCH --gres=$GRES
#SBATCH --output=$OUT_FILE
#SBATCH --error=${OUT_DIR}/%j.out # standard error output file named after job ID
#SBATCH --ntasks-per-node=$NTASKS
${MEM:+#SBATCH --mem=$MEM}

# PRINT JOB PARAMETERS TO ERROR OUTPUT FILE
(>&2 printf '%-17s: %s\n' 'SUBMIT DIRECTORY' '$(pwd)')
(>&2 printf '%-17s: %s\n' 'JOB NAME' '$JOB_NAME')
(>&2 printf '%-17s: %s\n' 'PARTITION' '$PARTITION')
(>&2 printf '%-17s: %s\n' 'GRES' '$GRES')
(>&2 printf '%-17s: %s\n' 'OUTPUT FILE' '$OUT_FILE')
(>&2 printf '%-17s: %s\n' 'TASKS PER NODE' '$NTASKS')
(>&2 printf '%-17s: %s\n' 'MEMORY REQUIRED' '$MEM')
(>&2 echo)
(>&2 printf '%-17s: %s\n' 'SCRIPT' '$SCRIPT')
(>&2 echo)

$SCRIPT # launching script
EOF
