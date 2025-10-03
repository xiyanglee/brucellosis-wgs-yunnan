#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /home/apps/anaconda3/envs/panaroo

set -uex

source _data_dir.sh

# nohup bash 17_3_panaroo_strict_threshold.sh >> _log/panaroo_strict_threshold.log 2>&1 &


panaroo \
    -i ${DATA_DIR}/_all_gff/*.gff \
    -o ${PANAROO_75_OUT_DIR} \
    --clean-mode strict \
    --family_threshold 0.75 \
    -t 12

panaroo \
    -i ${DATA_DIR}/_all_gff/*.gff \
    -o ${PANAROO_80_OUT_DIR} \
    --clean-mode strict \
    --family_threshold 0.80 \
    -t 12

panaroo \
    -i ${DATA_DIR}/_all_gff/*.gff \
    -o ${PANAROO_85_OUT_DIR} \
    --clean-mode strict \
    --family_threshold 0.85 \
    -t 12

panaroo \
    -i ${DATA_DIR}/_all_gff/*.gff \
    -o ${PANAROO_90_OUT_DIR} \
    --clean-mode strict \
    --family_threshold 0.90 \
    -t 12

panaroo \
    -i ${DATA_DIR}/_all_gff/*.gff \
    -o ${PANAROO_95_OUT_DIR} \
    --clean-mode strict \
    --family_threshold 0.95 \
    -t 12


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"