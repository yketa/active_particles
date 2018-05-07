#! /bin/bash

if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

DATA_DIRECTORY="/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}"
export DATA_DIRECTORY=$DATA_DIRECTORY

TIME=${TIME-1}
export INITIAL_FRAME=${INITIAL_FRAME-`/home/yketa/miniconda3/bin/python3.6 <<EOF
print(int(($(/home/yketa/bin/_colmig_DPD_P_A_data $DATA N_steps)//$(/home/yketa/bin/_colmig_DPD_P_A_data $DATA period_dump))/2))
EOF
`}
export INTERVAL_MAXIMUM=${INTERVAL_MAXIMUM-1}
export N_CASES=${N_CASES-`/home/yketa/miniconda3/bin/python3.6 <<EOF
from numpy import sqrt
N = $(/home/yketa/bin/_colmig_DPD_P_A_data $DATA N)
print(int(sqrt(N)) + (1 - int(sqrt(N))%2))
EOF
`}
export R_CUT=${R_CUT-2}
export SIGMA=${SIGMA-2}

export ENDPOINT=${ENDPOINT-False}

export SHOW='False'

TIME_LIST=( $TIME )
. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
if (( ${#TIME_LIST[@]} == 1 )); then
	CSS_TIME=Css_I$(float_to_letters $INITIAL_FRAME)_T$(float_to_letters $TIME)_M$(float_to_letters $INTERVAL_MAXIMUM)_C$(float_to_letters $N_CASES)_RCUT$(float_to_letters $R_CUT)_SIGM$(float_to_letters $SIGMA)
else
	CSS_TIME=Css_I$(float_to_letters $INITIAL_FRAME)_M$(float_to_letters $INTERVAL_MAXIMUM)_C$(float_to_letters $N_CASES)_RCUT$(float_to_letters $R_CUT)_SIGM$(float_to_letters $SIGMA)
fi

mkdir -p ${DATA_DIRECTORY}/out
output_file=${DATA_DIRECTORY}/out/${CSS_TIME}.out
echo "TIME: $TIME" > $output_file

sbatch --job-name=${CSS_TIME}_${DATA} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#!/bin/bash
$(OUTPUT=/home/yketa/hoomd/colmig_DPD_P_A/sub/out/${CSS_TIME}_${DATA}.%j.out /home/yketa/hoomd/colmig_DPD_P_A/sub/launch/header.sh)

for TIME in $TIME; do
	export TIME=\$TIME
	/home/yketa/bin/_colmig_DPD_P_A_Css >> $output_file
done
EOF
