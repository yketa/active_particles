if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

DATA_DIRECTORY="/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}"
export DATA_DIRECTORY=$DATA_DIRECTORY

export INITIAL_FRAME=${INITIAL_FRAME-`/home/yketa/miniconda3/bin/python3.6 <<EOF
print(int(($(/home/yketa/bin/_colmig_DPD_P_A_data $DATA N_steps)//$(/home/yketa/bin/_colmig_DPD_P_A_data $DATA period_dump))/2))
EOF
`}
export INTERVAL_MAXIMUM=${INTERVAL_MAXIMUM-1}
export N_CASES=${N_CASES-`/home/yketa/miniconda3/bin/python3.6 <<EOF
from numpy import sqrt
N = $(/home/yketa/bin/_colmig_DPD_P_A_data $DATA N)
print(int(np.sqrt(N)) + (1 - int(np.sqrt(N))%2))
EOF
`}

export SHOW='False'
export SAVE='True'

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
VARN_TIME=varN_I$(float_to_letters $INITIAL_FRAME)_M$(float_to_letters $INTERVAL_MAXIMUM)_C$(float_to_letters $N_CASES)

mkdir -p ${DATA_DIRECTORY}/out
output_file=${DATA_DIRECTORY}/out/${VARN_TIME}.out
> $output_file

sbatch --job-name=${VARN_TIME}_${DATA} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/${VARN_TIME}_${DATA}.%j.out
#SBATCH --ntasks-per-node 1

/home/yketa/miniconda3/bin/python3.6 /home/yketa/hoomd/scripts/varN.py >> $output_file
EOF

