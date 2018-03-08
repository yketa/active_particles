if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

export DATA_DIRECTORY=/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}
export PARAMETERS_FILE=${PARAMETERS_FILE-${DATA_DIRECTORY}/param.pickle}
export POSITION_FILE=${POSITION_FILE-${DATA_DIRECTORY}/position.csv}

export INITIAL_FRAME=${INITIAL_FRAME-0}
export SNAP_MAXIMUM=${SNAP_MAXIMUM-10}
export SNAP_PERIOD=${SNAP_PERIOD-1}

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
MSD_PAR=msd_I$(float_to_letters $INITIAL_FRAME)_M$(float_to_letters $SNAP_MAXIMUM)_P$(float_to_letters $SNAP_PERIOD)

mkdir -p ${DATA_DIRECTORY}/out
output_file=${DATA_DIRECTORY}/out/${MSD_PAR}.out
> $output_file

sbatch --job-name=${MSD_PAR}_${DATA} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/${MSD_PAR}_${DATA}.%j.out
#SBATCH --ntasks-per-node 1

/home/yketa/miniconda3/bin/python3.6 /home/yketa/hoomd/scripts/msd.py >> $output_file
EOF

