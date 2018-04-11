#! /bin/bash

if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

export DATA_DIRECTORY=/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}

export FRAME=${FRAME-`/home/yketa/miniconda3/bin/python3.6 <<EOF
print(int(($(/home/yketa/bin/_colmig_DPD_P_A_data $DATA N_steps)//$(/home/yketa/bin/_colmig_DPD_P_A_data $DATA period_dump))/2))
EOF
`}
export DT=${DT--1}
export MAX_BOX_SIZE=${MAX_BOX_SIZE-$(/home/yketa/bin/_colmig_DPD_P_A_data $DATA box_size)}

export MOVE=${MOVE-False}

export SHOW=False
export SAVE=True

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
FRA_PAR=frame_F$(float_to_letters $FRAME)_T$(float_to_letters $(($DT>0?$DT:0)))_RMAX$(float_to_letters $MAX_BOX_SIZE)_MOVE$(float_to_letters $MOVE)

mkdir -p /home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}/out
output_file=/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}/out/${FRA_PAR}.out
> $output_file

sbatch --job-name=${FRA_PAR}_${DATA} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/${FRA_PAR}_${DATA}.%j.out
#SBATCH --ntasks-per-node 1

/home/yketa/miniconda3/bin/python3.6 /home/yketa/bin/_colmig_DPD_P_A_frame >> $output_file
EOF

