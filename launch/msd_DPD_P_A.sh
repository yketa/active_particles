if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

export DATA_DIRECTORY=/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}

output_file=${DATA_DIRECTORY}/msd.out
> $output_file

sbatch --job-name=msd_${DATA} <<EOF
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --output /home/yketa/hoomd/colmig_DPD_P_A/sub/out/msd_${DATA}.%j.out
#SBATCH --ntasks-per-node 1

/home/yketa/miniconda3/bin/python3.6 /home/yketa/hoomd/scripts/msd.py >> $output_file
EOF

