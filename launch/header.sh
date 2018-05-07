#! /bin/bash

PARTITION=${PARTITION-gpu}
GRES=${GRES-gpu:k80:1}
OUTPUT=${OUTPUT-${PARTITION}.out}
NTASKS=${NTASKS-1}

cat <<EOF
#SBATCH --partition=$PARTITION
#SBATCH --gres=$GRES
#SBATCH --output $OUTPUT
#SBATCH --ntasks-per-node $NTASKS
EOF

