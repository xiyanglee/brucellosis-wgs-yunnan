#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /public/apps/roary-3.13.0

set -uex

source _data_dir.sh

# nohup bash 5_roary.sh >> _log/roary.log 2>&1 &


## GFF files must have unique basenames.
## $(ls ${PROKKA_OUT_DIR} | xargs -I {} echo ${PROKKA_OUT_DIR}/{}/prokka_prediction.gff)
ls ${PROKKA_OUT_DIR} | while read FILE_PREFIX
do
    cp ${PROKKA_OUT_DIR}/${FILE_PREFIX}/prokka_prediction.gff ${DATA_DIR}/_all_gff/${FILE_PREFIX}.gff
done

time roary \
    -i 95 \
    -cd 100 \
    -e \
    -n \
    -p 20 \
    -r \
    -v \
    -f ${ROARY_OUT_DIR} \
    ${DATA_DIR}/_all_gff/*.gff


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"