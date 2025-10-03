#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

set -uex

source _data_dir.sh

module load fastp/0.23.4

# nohup bash 3_fastp.sh >> _log/fastp.log 2>&1 &

cat ${FILE_LIST} | awk 'NR>=2{print $0}' | while read LINE
do
    SAMPLE_NAME=$(echo ${LINE} | cut -d ' ' -f 1)
    R1_PATH=$(echo ${LINE} | cut -d ' ' -f 2)
    R2_PATH=$(echo ${LINE} | cut -d ' ' -f 3)
    fastp \
        -i ${R1_PATH} \
        -I ${R2_PATH} \
        --html ${FASTP_OUT_DIR}/${SAMPLE_NAME}_cleaned.html \
        --json ${FASTP_OUT_DIR}/${SAMPLE_NAME}_cleaned.json \
        -o ${FASTP_OUT_DIR}/${SAMPLE_NAME}_clean_1.fq.gz \
        -O ${FASTP_OUT_DIR}/${SAMPLE_NAME}_clean_2.fq.gz \
        --low_complexity_filter \
        --trim_poly_g \
        --trim_poly_x \
        --thread 16 \
        1>${FASTP_OUT_DIR}/${SAMPLE_NAME}.log 2>&1
done


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"