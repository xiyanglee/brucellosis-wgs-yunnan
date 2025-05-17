#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

# -u: 检查脚本内变量，若有未定义变量则停止运行
# -e: 若脚本遇到错误则停止执行
# -x: 显示脚本执行过程
set -uex

source _data_dir.sh

module load fastqc/0.12.1

# nohup bash 1_fastqc.sh >> _log/fastqc.log 2>&1 &

fastqc \
    --noextract \
    --outdir ${FASTQC_OUT_DIR} \
    --threads 16 \
    $(cat ${FILE_LIST} | awk -v OFS='\n' 'NR>=2{print $2,$3}')

set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"