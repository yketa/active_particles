if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

DATA_DIRECTORY="/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}"
export VELOCITY_FILE="${DATA_DIRECTORY}/velocity.csv"
export PARAMETERS_FILE="${DATA_DIRECTORY}/param.pickle"
export POSITION_FILE="${DATA_DIRECTORY}/position.csv"

mkdir -p ${DATA_DIRECTORY}/out
output_file=${DATA_DIRECTORY}/out/pos.out
> $output_file

sbatch --job-name=pos_${DATA} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/pos_${DATA}.%j.out
#SBATCH --ntasks-per-node 1

/home/yketa/bin/_pos >> $output_file
EOF

