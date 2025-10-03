#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /home/apps/anaconda3/envs/panaroo

set -uex

source _data_dir.sh

# nohup bash 17_1_panaroo_sensitive_mode.sh >> _log/panaroo_sensitive.log 2>&1 &


panaroo \
    -i ${DATA_DIR}/_all_gff/*.gff \
    -o ${PANAROO_SENSITIVE_OUT_DIR} \
    --clean-mode sensitive \
    -t 12


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"