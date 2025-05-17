#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

set +x
source /home/apps/anaconda3/bin/activate
conda activate /public/apps/multiqc-1.23
set -uex

source _data_dir.sh

# nohup bash 2_multiqc.sh >> _log/multiqc.log 2>&1 &

multiqc \
    --outdir ${MULTIQC_OUT_DIR} \
    ${FASTQC_OUT_DIR}

set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"