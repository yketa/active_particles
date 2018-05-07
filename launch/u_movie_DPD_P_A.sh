#! /bin/bash

if [[ "$DATA" == "" ]]; then
	echo "No data name submitted."
	exit 0
fi

export DATA=$DATA

export INITIAL_FRAME=${INITIAL_FRAME-0}
export FINAL_FRAME=${FINAL_FRAME-10000}
export FRAME_PERIOD=${FRAME_PERIOD-1}
export FRAME_MAXIMUM=${FRAME_MAXIMUM-1000}
export TIME_SCALE=${TIME_SCALE-LOG}

export V_MIN=${V_MIN--1}
export V_MAX=${V_MAX--1}

export MAX_BOX_SIZE=${MAX_BOX_SIZE--1}

export MOVE=${MOVE-False}

. /home/yketa/exponents.sh # exporting letters expressions and float conversion functions
MOV_PAR=u_mov_I$(float_to_letters $INITIAL_FRAME)_M$(float_to_letters $FRAME_MAXIMUM)_P$(float_to_letters $FRAME_PERIOD)_MOVE$(float_to_letters $MOVE)

MOV=u_${DATA}_I$(float_to_letters $INITIAL_FRAME)_M$(float_to_letters $FRAME_MAXIMUM)_P$(float_to_letters $FRAME_PERIOD)_MOVE$(float_to_letters $MOVE)
export MOVIE_DIRECTORY=${MOVIE_DIRECTORY-/home/yketa/hoomd/colmig_DPD_P_A/movie/${MOV}}

mkdir -p /home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}/out
output_file=/home/yketa/hoomd/colmig_DPD_P_A/data/${DATA}/out/${MOV_PAR}.out
> $output_file

sbatch --job-name=${MOV_PAR}_${DATA} <<EOF
#!/bin/bash
$(OUTPUT=/home/yketa/hoomd/colmig_DPD_P_A/sub/out/${MOV_PAR}_${DATA}.%j.out /home/yketa/hoomd/colmig_DPD_P_A/sub/launch/header.sh)

/bin/bash /home/yketa/bin/_colmig_DPD_P_A_u_makemovie >> $output_file
EOF
