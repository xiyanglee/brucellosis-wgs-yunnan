#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

set -uex

source _data_dir.sh

# nohup bash 10_1_snippy_contig_alignment.sh >> _log/snippy_contig_alignment.log 2>&1 &


cat ${FILE_LIST} | awk 'NR>=2{print $1}' | while read SAMPLE_NAME
do
    cp ${UNICYCLER_OUT_DIR}/${SAMPLE_NAME}/assembly.fasta ${UNICYCLER_ALL_FNA}/${SAMPLE_NAME}-assembly.fasta
done

cat ${SOURCE_DIR}/assembly_genbank/assembly_filepath.tsv | awk 'NR>=2{print $0}' | while read LINE
do
    SAMPLE_NAME=$(echo ${LINE} | cut -d ' ' -f 1)
    FILE_PREFIX=$(echo ${LINE} | cut -d ' ' -f 2)
    cp ${SOURCE_DIR}/assembly_genbank/genomes/${FILE_PREFIX}/${FILE_PREFIX}_genomic.fna.gz ${UNICYCLER_ALL_FNA}/${SAMPLE_NAME}-genebank.fna.gz
done

# [08:43:34] Shredding /data/data/unicycler_and_source_fna/GCA_000022625.1-genebank.fna.gz into pseudo-reads

# ------------- EXCEPTION -------------
# MSG: The sequence does not appear to be FASTA format (lacks a descriptor line '>')
# STACK Bio::SeqIO::fasta::next_seq /usr/share/perl5/Bio/SeqIO/fasta.pm:136
# STACK toplevel /snippy-4.6.0/bin/snippy:255
# -------------------------------------

for file in *.gz; do [ -f "$file" ] && gzip -d "$file"; done

ls ${UNICYCLER_ALL_FNA}/* | grep 'GCA_' | grep -v 'GCA_000007125' | while read INPUT_FILE
do
    FILE_NAME=$(basename ${INPUT_FILE} .txt)
    SAMPLE_NAME=$(echo ${FILE_NAME} | cut -d '-' -f 1)
    docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) snippy:4.6.0_SC2 \
        snippy \
            --outdir /data/data/snippy_contig_alignment/${SAMPLE_NAME} \
            --ctgs /data/data/unicycler_and_source_fna/${FILE_NAME} \
            --ref /data/source/reference_genome/GCF_000007125.1.gbff \
            --cpus 12
done

docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) snippy:4.6.0_SC2 \
    snippy \
        --outdir /data/data/snippy_contig_alignment/GCA_000369945.1 \
        --ctgs /data/source/out_group/GCA_000369945.1.fna \
        --ref /data/source/reference_genome/GCF_000007125.1.gbff \
        --cpus 12

docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) -w /data/data/snippy_contig_alignment snippy:4.6.0_SC2 \
    snippy-core --ref /data/source/reference_genome/GCF_000007125.1.gbff $(ls -1 ${SNIPPY_CONTIG_OUT})


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"