#! /bin/bash

PARTITION=${PARTITION-gpu}
GPU=${GPU-gpu:k80:1}
OUTPUT=${OUTPUT-${PARTITION}.out}
NTASKS=${NTASKS-1}

cat <<EOF
#SBATCH --partition=$PARTITION
#SBATCH --gres=$GPU
#SBATCH --output $OUTPUT
#SBATCH --ntasks-per-node $NTASKS
EOF

