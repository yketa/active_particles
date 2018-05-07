#! /bin/bash

export N=${NUMBER-2000}
export N_STEPS=${N_STEPS-1e4}
export PERIOD_DUMP=${PERIOD_DUMP-1e2}
export TIME_STEP=${TIME_STEP-1e-2}

export DENSITY=$DENSITY
export SELF_PROPULSION_SPEED=$SELF_PROPULSION_SPEED
export ROTATION_DIFFUSION=$ROTATION_DIFFUSION

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
DATA=D$(float_to_letters $DENSITY)_V$(float_to_letters $SELF_PROPULSION_SPEED)_R$(float_to_letters $ROTATION_DIFFUSION)_N$(float_to_letters $N) # name of the data
sim_name=${DATA}_L$(float_to_letters ${OVERWRITE-$(ls -d /home/yketa/hoomd/colmig_DPD_P_A/data/*/ | grep $DATA | wc -l)})

export DATA_DIRECTORY=/home/yketa/hoomd/colmig_DPD_P_A/data/${sim_name}

output_file=/home/yketa/hoomd/colmig_DPD_P_A/out/${sim_name}.out
> $output_file

sbatch --job-name=$sim_name ${CHAIN:+-d afterok:$CHAIN} <<EOF
#!/bin/bash
$(OUTPUT=/home/yketa/hoomd/colmig_DPD_P_A/sub/out/${sim_name}.%j.out /home/yketa/hoomd/colmig_DPD_P_A/sub/launch/header.sh)

/home/yketa/miniconda3/bin/python3.6 /home/yketa/hoomd/scripts/collective_migration_DPD_polydisperse_active.py >> $output_file
EOF
