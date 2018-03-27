export N=${NUMBER-2000}
export N_STEPS=${N_STEPS-1e4}
export PERIOD_DUMP=${PERIOD_DUMP-1e2}
export TIME_STEP=${TIME_STEP-1e-2}
export SPRING_CONSTANT=${SPRING_CONSTANT-1}

export DENSITY=$DENSITY
export TEMPERATURE=$TEMPERATURE

export SELF_PROPULSION_SPEED=0
export ROTATION_DIFFUSION=0

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
DATA=D$(float_to_letters $DENSITY)_K$(float_to_letters $TEMPERATURE)_S$(float_to_letters $SPRING_CONSTANT)_N$(float_to_letters $N) # name of the data
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

/home/yketa/miniconda3/bin/python3.6 /home/yketa/hoomd/scripts/collective_migration_DPD_polydisperse_active.py >> $output_file
EOF

