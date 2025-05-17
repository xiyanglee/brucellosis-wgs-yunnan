#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
# diamond-0.9.14

set -uex

source _data_dir.sh

# nohup bash 11_1_diamond_VFDB.sh >> _log/diamond_VFDB.log 2>&1 &



# diamond makedb --in ${SOURCE_DIR}/VFDB/VFDB_setB_pro.fas.gz --db ${DIAMOND_DIR}/VFDB_setB

diamond blastx \
    --db ${DIAMOND_DIR}/VFDB_setB \
    --query ${ROARY_OUT_DIR}/pan_genome_reference.fa \
    -o ${DIAMOND_DIR}/diamond_SetB_out.tsv \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --min-score 60 \
    --block-size 40.0 \
    --index-chunks 1

# cat /public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly/data/roary_result/pan_genome_reference.fa | grep '^>' | sed 's/^>//g' > /public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly/data/diamond/pan_genome_reference_id.txt

set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"