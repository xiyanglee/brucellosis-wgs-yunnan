#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

# -e：exit on error —— 遇到错误立即退出
# -u：treat unset variables as error —— 使用未定义变量时报错
# -o pipefail：fail if any command in a pipeline fails
# 在 A | B | C 这样的管道中，如果 A 或 B 失败（返回非 0），整体会被视为失败
set -euo pipefail

source _data_dir.sh

# nohup bash 15_1_fastANI.sh >> _log/fastANI.log 2>&1 &


THREADS=14
FNA_LIST=${FASTANI_OUT}/fna_list.txt
PAIR_LIST=${FASTANI_OUT}/pair_list.tsv
RESULT_DIR=${FASTANI_OUT}/results
RESULT_FILE=${FASTANI_OUT}/fastani_all_results.out

# 创建输出目录
mkdir -p "${RESULT_DIR}"

# 获取所有fna文件
ls "${PROKKA_OUT_DIR}" | while read FILE_PREFIX; do
    echo -e "${FILE_PREFIX}\t${PROKKA_OUT_DIR}/${FILE_PREFIX}/prokka_prediction.fna"
done > "${FNA_LIST}"

# 生成所有非冗余组合（i < j）
awk 'NR==FNR{a[NR]=$1; p[NR]=$2; next}
     {for(i=1;i<FNR;i++) print a[i] "\t" p[i] "\t" $1 "\t" $2}' "${FNA_LIST}" "${FNA_LIST}" > "${PAIR_LIST}"

# 并发执行 fastANI，每个输出写入独立文件
cat "${PAIR_LIST}" | \
xargs -n 4 -P ${THREADS} -I{} bash -c '
    LINE="{}"
    PREFIX_A=$(echo "$LINE" | cut -f1)
    FILE_A=$(echo "$LINE" | cut -f2)
    PREFIX_B=$(echo "$LINE" | cut -f3)
    FILE_B=$(echo "$LINE" | cut -f4)
    OUTFILE="'${RESULT_DIR}'/${PREFIX_A}_vs_${PREFIX_B}.out"

    if [[ ! -s "$OUTFILE" ]]; then
        /home/xiyangli/alto/software/fastANI -q "$FILE_A" -r "$FILE_B" -o "$OUTFILE"
    fi
'

# 合并所有结果文件
cat ${RESULT_DIR}/*.out > "${RESULT_FILE}"


set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"