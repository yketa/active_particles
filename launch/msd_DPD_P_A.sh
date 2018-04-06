if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

export DATA_DIRECTORY=/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}
export PARAMETERS_FILE=${PARAMETERS_FILE-${DATA_DIRECTORY}/param.pickle}
export UNWRAPPED_FILE=${UNWRAPPED_FILE-${DATA_DIRECTORY}/trajectory.dat}

INITIAL_FRAME=${INITIAL_FRAME-`/home/yketa/miniconda3/bin/python3.6 <<EOF
print(int(($(/home/yketa/bin/_colmig_DPD_P_A_data $DATA N_steps)//$(/home/yketa/bin/_colmig_DPD_P_A_data $DATA period_dump))/2))
EOF
`}
export SNAP_MAXIMUM=${SNAP_MAXIMUM-1}
export SNAP_PERIOD=${SNAP_PERIOD-1}

INITIAL_FRAME_LIST=( $INITIAL_FRAME )
. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
if (( ${#INITIAL_FRAME_LIST[@]} == 1 )); then
	MSD_PAR=msd_I$(float_to_letters $INITIAL_FRAME)_M$(float_to_letters $SNAP_MAXIMUM)_P$(float_to_letters $SNAP_PERIOD)
else
	MSD_PAR=msd_M$(float_to_letters $SNAP_MAXIMUM)_P$(float_to_letters $SNAP_PERIOD)
fi

mkdir -p ${DATA_DIRECTORY}/out
output_file=${DATA_DIRECTORY}/out/${MSD_PAR}.out
echo "INITIAL_FRAME: $INITIAL_FRAME"> $output_file

sbatch --job-name=${MSD_PAR}_${DATA} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/${MSD_PAR}_${DATA}.%j.out
#SBATCH --ntasks-per-node 1

for INITIAL_FRAME in $INITIAL_FRAME; do
	export INITIAL_FRAME=\$INITIAL_FRAME
	/home/yketa/miniconda3/bin/python3.6 /home/yketa/hoomd/scripts/msd.py >> $output_file
done
EOF

