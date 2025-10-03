#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

source /home/apps/anaconda3/bin/activate
conda activate eggnog-mapper

set -uex

source _data_dir.sh

# nohup bash 17_5_emapper_diamond_panaroo_strict_pan.sh >> _log/emapper_diamond_panaroo_strict_pan.log 2>&1 &

## 其他参数
# --mp_start_method [fork,spawn,forkserver] 控制 Python multiprocessing 启动方法。仅当 multiprocessing 未在您的作系统中正常运行时，才使用

## 比对参数
# --pident FLOAT         仅报告等于或高于此身份阈值百分比 (0.0-100.0) 的对齐；如果 -m hmmer 则无效
# --evalue FLOAT         仅报告等于或高于此 E 值阈值的对齐。默认 0.001
# --score FLOAT          仅报告等于或高于此 bit score 阈值的对齐。默认 None
# --query_cover FLOAT    仅报告等于或高于此查询覆盖率分数阈值 (0.0 - 100.0) 的对齐。默认 None。如果 -m hmmer 则无效
# --subject_cover FLOAT  仅报告等于或高于此目标（eggNOG 序列）覆盖率分数阈值 (0.0 - 100.0) 的比对。默认 None。如果 -m hmmer 则无效

## 基本设置
# --override             覆盖输出文件（如果存在）。默认情况下，如果检测到冲突文件，则会中止执行。
# --dmnd_ignore_warnings Diamond 的 --ignore-warnings 选项。使用它来避免 Diamond 因警告而停止执行
# --itype CDS
# --translate            仅在 CDS 模式，先翻译为蛋白再比对
# --tax_scope            限制注释范围
# --target_orthologs all 所有类型的直系/旁系同源都考虑
# --go_evidence all      获取全部 GO 注释，包括电子推导
# --pfam_realign none    不做结构域重新比对 (节省时间)
# --report_orthologs     生成 .orthologs 文件
# --excel                额外导出 Excel 格式结果

    # --score 50 \
    # --pident 30 \
    # --query_cover 20 \
    # --subject_cover 20 \
emapper.py \
    -i ${PANAROO_PAN_OUT_DIR}/pan_genome_reference.fa \
    -o pan_genome_annotation \
    --output_dir ${PANAROO_EMAPPER_DIAMOND_PAN_OUT_DIR} \
    --data_dir /public/data_2/lxy/alto/db/eggnog_mapper_bacteria/ \
    --cpu 14 \
    --override \
    -m diamond \
    --dmnd_ignore_warnings \
    --evalue 0.001 \
    --itype CDS \
    --translate \
    --tax_scope 2 \
    --target_orthologs all \
    --go_evidence all \
    --pfam_realign none \
    --report_orthologs \
    --excel



set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"