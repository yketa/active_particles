if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

DATA_DIRECTORY="/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}"
export DATA_DIRECTORY=$DATA_DIRECTORY
export TIME=${TIME-1}

export SHOW='False'

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
CUU_TIME=Cuu_T$(float_to_letters $TIME)

mkdir -p ${DATA_DIRECTORY}/out
output_file=${DATA_DIRECTORY}/out/${CUU_TIME}.out
> $output_file

sbatch --job-name=${CUU_TIME}_${DATA} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/${CUU_TIME}_${DATA}.%j.out
#SBATCH --ntasks-per-node 1

/home/yketa/bin/_colmig_DPD_P_A_Cuu >> $output_file
EOF

