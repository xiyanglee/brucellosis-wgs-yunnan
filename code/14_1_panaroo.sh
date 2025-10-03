#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /home/apps/anaconda3/envs/panaroo

set -uex

source _data_dir.sh

# nohup bash 14_1_panaroo.sh >> _log/panaroo.log 2>&1 &


# ls ${PROKKA_OUT_DIR} | while read FILE_PREFIX
# do
#     cp ${PROKKA_OUT_DIR}/${FILE_PREFIX}/prokka_prediction.gff ${DATA_DIR}/_all_gff/${FILE_PREFIX}.gff
# done
panaroo \
    -i ${DATA_DIR}/_all_gff/*.gff \
    -o ${PANAROO_OUT_DIR} \
    --clean-mode strict \
    -a core \
    -t 12


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"