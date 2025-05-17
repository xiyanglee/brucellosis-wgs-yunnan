#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate /public/apps/raxml-8.2.13

set -uex

source _data_dir.sh

# nohup bash 8_3_raxml.sh >> _log/raxml_8_3.log 2>&1 &


cd ${REMOVED_GAP_DIR}
# -s 用于建树的序列
# -n 输出的树文件的后缀名
# -m 模型设置
# -T 线程数
# -N 自举次数，一般论文要求 1000
# -p (小写) 一个用于简约推断的随机数，随机数可以自己设 (如 12345)
# -f a 快速自举分析
# -x (小写) 指定一个整数 (随机数，可以随便取比如 12345) 并启用快速 bootstrapping
# -k (小写) 指定 bootstrapped 树应该输出分支长度, 默认关闭
raxmlHPC-PTHREADS-AVX2 \
	-s ${REMOVED_GAP_DIR}/core_gene_snp_filtered_unique.phy \
	-n raxml_out \
	-m GTRGAMMAI \
	-T 12 \
	-N 1000 \
	-p 12345 \
	-x 12345 \
	-f a 


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"