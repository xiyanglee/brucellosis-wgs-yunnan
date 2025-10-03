#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /public/apps/emboss-6.6.0

set -uex

source _data_dir.sh

module load seqtk/1.4

# nohup bash 7_0_seqtk.sh >> _log/seqtk.log 2>&1 &


cat ${FILE_LIST} | awk 'NR>=2{print $0}' | cut -f 1 | while read SAMPLE_NAME
do
	echo ${SAMPLE_NAME} >> ${TMP_DIR}/sample_list_yunnan.txt
done

time seqtk subseq \
	${ROARY_OUT_DIR}/core_gene_alignment.aln \
	${TMP_DIR}/sample_list_yunnan.txt \
	> ${TMP_DIR}/core_gene_alignment_yunnan.aln

time seqret \
	-sequence ${TMP_DIR}/core_gene_alignment_yunnan.aln \
	-outseq ${TMP_DIR}/core_gene_alignment_yunnan.nex \
	-osformat nexus


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"