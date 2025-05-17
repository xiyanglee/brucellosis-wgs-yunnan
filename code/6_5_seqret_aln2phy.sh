#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /public/apps/emboss-6.6.0

set -uex

source _data_dir.sh

# nohup bash 6_5_seqret_aln2phy.sh >> _log/seqret_aln2phy.log 2>&1 &


ls ${PROKKA_OUT_DIR} | grep GCA_ | awk -v OFS=";" '{print $1,$1}' \
	| sed 's/;GCA_/;/g' | sed 's/\.[0-9]$//g' | sed 's/;/\t/g' \
	> ${TRIMAL_OUT_DIR}/phylip_sample_name_genebank.tsv

ls ${PROKKA_OUT_DIR} | grep GCA_ | awk -v OFS=";" '{print $1,$1}' \
	| sed 's/;GCA_/;/' | sed 's/\.[0-9]$//' | sed 's/\./\\\./' | sed 's/^/s\//' | sed 's/;/\/>/g' | sed 's/GCA_/\^>GCA_/' | sed 's/$/\//' \
	> ${TRIMAL_OUT_DIR}/phylip_sample_name_genebank.sed

sed -i.bak -f ${TRIMAL_OUT_DIR}/phylip_sample_name_genebank.sed ${TRIMAL_OUT_DIR}/core_gene_alignment_trim_unique.aln

time seqret \
	-sequence ${TRIMAL_OUT_DIR}/core_gene_alignment_trim_unique.aln \
	-outseq ${TRIMAL_OUT_DIR}/core_gene_alignment_trim_unique.phy \
	-osformat phylip


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"