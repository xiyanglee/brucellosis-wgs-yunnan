#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /public/apps/prokka-1.14.6

set -uex

source _data_dir.sh

# nohup bash 4_prokka.sh >> _log/prokka.log 2>&1 &


cat ${FILE_LIST} | awk 'NR>=2{print $1}' | while read SAMPLE_NAME
do
    ## use unicycler results
    cp ${UNICYCLER_OUT_DIR}/${SAMPLE_NAME}/assembly.fasta ${DATA_DIR}/_genebank_fasta/${SAMPLE_NAME}-assembly.fasta
done

cat ${SOURCE_DIR}/assembly_genbank/assembly_filepath.tsv | awk 'NR>=2{print $0}' | while read LINE
do
    SAMPLE_NAME=$(echo ${LINE} | cut -d ' ' -f 1)
    FILE_PREFIX=$(echo ${LINE} | cut -d ' ' -f 2)
    gunzip -c ${SOURCE_DIR}/assembly_genbank/genomes/${FILE_PREFIX}/${FILE_PREFIX}_genomic.fna.gz > ${DATA_DIR}/_genebank_fasta/${SAMPLE_NAME}-genebank.fasta
done

## https://github.com/tseemann/prokka
# ls ${MLST_OUT_DIR}/unicycler_fasta/* | while read INPUT_FILE_PATH
## prokka 似乎无法读取 gz 压缩文件
ls ${DATA_DIR}/_genebank_fasta/* | while read INPUT_FILE_PATH
do
    FILE_PREFIX=$(echo ${INPUT_FILE_PATH##*/} | cut -d '-' -f 1)

    # --outdir 无需提前创建，否则将报错；使用已有路径需添加参数 --force
    time prokka \
        --kingdom Bacteria \
        --genus Brucella \
        --species melitensis \
        --locustag PROKKATAG-${FILE_PREFIX} \
        --cpus 20 \
        --mincontiglen 500 \
        --prefix prokka_prediction \
        --outdir ${PROKKA_OUT_DIR}/${FILE_PREFIX} \
        ${INPUT_FILE_PATH}
done

# rm -rf ${DATA_DIR}/_genebank_fasta/


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"