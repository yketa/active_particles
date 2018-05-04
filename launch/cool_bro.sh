#! /bin/bash

export N=${NUMBER-1e5}
export N_STEPS=${N_STEPS-1e4}
export PERIOD_DUMP=${PERIOD_DUMP-1e2}
export TIME_STEP=${TIME_STEP-1e-2}

export DENSITY=${DENSITY-0.8}
export TEMPERATURE_INITIAL=${TEMPERATURE_INITIAL-1}
export TEMPERATURE_FINAL=${TEMPERATURE_FINAL-0}
export PREPARATION_TIME=${PREPARATION_TIME-10}
export COOLING_RATE=${COOLING_RATE-1e-2}

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
DATA=D$(float_to_letters $DENSITY)_K0$(float_to_letters $TEMPERATURE_INITIAL)_K1$(float_to_letters $TEMPERATURE_FINAL)_CR$(float_to_letters $COOLING_RATE)_N$(float_to_letters $N) # name of the data
sim_name=${DATA}_L$(float_to_letters ${OVERWRITE-$(ls -d /home/yketa/hoomd/colmig_DPD_P_A/data/*/ | grep $DATA | wc -l)})

export DATA_DIRECTORY=/home/yketa/hoomd/colmig_DPD_P_A/data/${sim_name}

output_file=/home/yketa/hoomd/colmig_DPD_P_A/out/${sim_name}.out
> $output_file

sbatch --job-name=$sim_name ${CHAIN:+-d afterok:$CHAIN} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/${sim_name}.%j.out
#SBATCH --ntasks-per-node 1

/home/yketa/miniconda3/bin/python3.6 /home/yketa/hoomd/scripts/cooling_brownian.py >> $output_file
EOF
