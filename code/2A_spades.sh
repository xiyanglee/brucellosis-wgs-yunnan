#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

set -uex

source _data_dir.sh

module load SPAdes/4.0.0

# nohup bash 4_spades.sh >> _log/spades.log 2>&1 &

# https://ablab.github.io/spades/datatypes.html#assembling-long-illumina-paired-reads
# for 150bp reads SPAdes uses k-mer sizes 21, 33, 55, 77
cat ${FILE_LIST} | awk 'NR>=2{print $1}' | while read SAMPLE_NAME
do
    time spades.py \
        --isolate \
        -o ${SPADES_OUT_DIR}/${SAMPLE_NAME} \
        -t 20 \
        -m 120 \
        --tmp-dir /public/tmp \
        -1 ${FASTP_OUT_DIR}/${SAMPLE_NAME}_clean_1.fq.gz \
        -2 ${FASTP_OUT_DIR}/${SAMPLE_NAME}_clean_2.fq.gz \
        1>${SPADES_OUT_DIR}/${SAMPLE_NAME}.log \
        2>${SPADES_OUT_DIR}/${SAMPLE_NAME}.err
done


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"