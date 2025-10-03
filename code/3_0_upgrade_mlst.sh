#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /public/apps/mlst-2.23.0

set -uex

source _data_dir.sh

# nohup bash 5_mlst.sh >> _log/mlst.log 2>&1 &


ASSEMBLER=unicycler

FASTA_INPUT_DIR=${MLST_OUT_DIR}/${ASSEMBLER}_fasta
mkdir -p ${FASTA_INPUT_DIR}

cat ${FILE_LIST} | awk 'NR>=2{print $1}' | while read SAMPLE_NAME
do
    ## use unicycler results
    ln -s ${UNICYCLER_OUT_DIR}/${SAMPLE_NAME}/assembly.fasta ${FASTA_INPUT_DIR}/${SAMPLE_NAME}-genomic.fasta
done

cat ${SOURCE_DIR}/assembly_genbank/assembly_filepath.tsv | awk 'NR>=2{print $0}' | while read LINE
do
    SAMPLE_NAME=$(echo ${LINE} | cut -d ' ' -f 1)
    FILE_PREFIX=$(echo ${LINE} | cut -d ' ' -f 2)
    ln -s ${SOURCE_DIR}/assembly_genbank/genomes/${FILE_PREFIX}/${FILE_PREFIX}_genomic.fna.gz ${FASTA_INPUT_DIR}/${SAMPLE_NAME}-genomic.fna.gz
done

time mlst \
    --legacy \
    --scheme brucella \
    --nopath ${FASTA_INPUT_DIR}/* \
    --threads 16 \
    > ${MLST_OUT_DIR}/${ASSEMBLER}_mlst_st.tsv


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"