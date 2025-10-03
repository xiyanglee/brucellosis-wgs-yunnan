#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

set -uex

source _data_dir.sh

# nohup bash 13_2_snippy_coregenes_alignment.sh >> _log/snippy_coregenes_alignment.log 2>&1 &


# --ref /data/source/reference_genome/GCF_000007125.1.gbff \
ls ${COREGENE_MERGE_DIR}/* | grep -v 'GCA_000007125' | while read INPUT_FILE
do
    SAMPLE_NAME=$(basename ${INPUT_FILE} .fasta)
    FILE_NAME=${SAMPLE_NAME}.fasta
    docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) snippy:4.6.0_SC2 \
        snippy \
            --outdir /data/data/coregene_merge_snippy_alignment/${SAMPLE_NAME} \
            --ctgs /data/data/core_genes_fasta/${FILE_NAME} \
            --ref /data/data/core_genes_fasta/GCA_000007125.1.fasta \
            --cpus 12
done

docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) snippy:4.6.0_SC2 \
    snippy \
        --outdir /data/data/coregene_merge_snippy_alignment/GCA_000369945.1 \
        --ctgs /data/source/out_group/GCA_000369945.1.fna \
        --ref /data/data/core_genes_fasta/GCA_000007125.1.fasta \
        --cpus 12

docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) -w /data/data/coregene_merge_snippy_alignment snippy:4.6.0_SC2 \
    snippy-core --ref /data/data/core_genes_fasta/GCA_000007125.1.fasta $(ls -1 ${SNIPPY_COREGENE_OUT_DIR})

# cp coregene_merge_snippy_alignment/core.* coregene_merge_snippy_alignment_merge/

set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"