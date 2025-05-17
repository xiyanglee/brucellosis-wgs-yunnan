#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

# set -uex

source _data_dir.sh

# nohup bash 0_download_ncbi_assembly.sh >> _log/download_ncbi_assembly.log 2>&1 &


# echo -e "#assembly_accession\tfile_prefix" >> ${SOURCE_DIR}/assembly_genbank/assembly_filepath.tsv
# cat ${SOURCE_DIR}/assembly_genbank/assembly_summary_genbank.txt | grep 'Brucella melitensis' | cut -f 1,20 | while read LINE
# do
#     SAMPLE_NAME=$(echo ${LINE} | cut -d ' ' -f 1)
#     FTP_PATH=$(echo ${LINE} | cut -d ' ' -f 2)
#     DOWNLOAD_FILE_PREFIX=${FTP_PATH##*/}
#     echo -e "${SAMPLE_NAME}\t${DOWNLOAD_FILE_PREFIX}" >> ${SOURCE_DIR}/assembly_genbank/assembly_filepath.tsv
# done

cat ${SOURCE_DIR}/assembly_genbank/assembly_summary_genbank.txt | grep 'Brucella melitensis' | cut -f 1,20 | while read LINE
do
    SAMPLE_NAME=$(echo ${LINE} | cut -d ' ' -f 1)
    FTP_PATH=$(echo ${LINE} | cut -d ' ' -f 2)
    DOWNLOAD_FILE_PREFIX=${FTP_PATH##*/}

    mkdir -p ${SOURCE_DIR}/assembly_genbank/genomes/${DOWNLOAD_FILE_PREFIX}/
    # -c 表示断点可以续传
    # -t 重试次数 (default: 20)
    # --timeout 超时时间 s (default: 900)
    # wget -c ${FTP_PATH}/${DOWNLOAD_FILE} -O ${SOURCE_DIR}/assembly_genbank/genomes/${SAMPLE_NAME}_genomic.fna.gz
    # wget -c -t 100 ${FTP_PATH}/${DOWNLOAD_FILE_PREFIX}_genomic.fna.gz -P ${SOURCE_DIR}/assembly_genbank/genomes/${DOWNLOAD_FILE_PREFIX}/

    axel -n 16 -o ${SOURCE_DIR}/assembly_genbank/genomes/${DOWNLOAD_FILE_PREFIX}/ ${FTP_PATH}/${DOWNLOAD_FILE_PREFIX}_genomic.fna.gz
done

head -1 ${SOURCE_DIR}/assembly_genbank/assembly_summary_genbank.txt >> ${METADATA_DIR}/assembly_summary_genbank.tsv
cat ${SOURCE_DIR}/assembly_genbank/assembly_summary_genbank.txt | grep 'Brucella melitensis' >> ${METADATA_DIR}/assembly_summary_genbank.tsv


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"