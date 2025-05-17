#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /public/apps/unicycler-0.5.0

set -uex

source _data_dir.sh

module load blast/2.16.0
module load SPAdes/3.15.3
# module load SPAdes/4.0.0
# Error: Unicycler requires SPAdes v3.14.0 or later
# Dependencies:
#   Program       Version   Status      Path                 
#   spades.py     4.0.0     too new     /public/apps/SPAdes-4.0.0/bin/spades.py
# 
# finished abnormally, OS return value: -11
# Splitting kmer instances into 200 files using 20 threads. This might take a while.
# 
# module 'collections' has no attribute 'Hashable'
                  
# nohup bash 4_unicycler.sh >> _log/unicycler.log 2>&1 &

cat ${FILE_LIST} | awk 'NR>=2{print $1}' | while read SAMPLE_NAME
do
    ## https://github.com/rrwick/Unicycler
    ## some default args
    # --mode normal
    # --min_fasta_length 100
    # SPAdes assembly:
    # --min_kmer_frac 0.2
    # --max_kmer_frac 0.95
    # --kmer_count 8
    # --depth_filter 0.25
    time unicycler \
        -1 ${FASTP_OUT_DIR}/${SAMPLE_NAME}_clean_1.fq.gz \
        -2 ${FASTP_OUT_DIR}/${SAMPLE_NAME}_clean_2.fq.gz \
        -o ${UNICYCLER_OUT_DIR}/${SAMPLE_NAME} \
        -t 20 \
        --verbosity 2 \
        --keep 2 \
        --spades_options "-m 120 --tmp-dir /public/tmp" \
        1>${UNICYCLER_OUT_DIR}/${SAMPLE_NAME}.log \
        2>${UNICYCLER_OUT_DIR}/${SAMPLE_NAME}.err
done


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"