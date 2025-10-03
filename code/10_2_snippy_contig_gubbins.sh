#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

set -uex

source _data_dir.sh

# nohup bash 10_2_snippy_contig_gubbins.sh >> _log/snippy_contig_gubbins.log 2>&1 &


docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) snippy:4.6.0_SC2 \
    bash -c "snippy-clean_full_aln /data/data/snippy_contig_alignment/core.full.aln > /data/data/snippy_contig_gubbins/clean.full.aln"

# https://github.com/nickjcroucher/gubbins
# sudo apt-get install gubbins
# https://www.linuxfordevices.com/tutorials/linux/install-pthread-library-linux
# sudo apt-get install libpthread-stubs0-dev
cd ${SUB_PROJECT_DIR}/data/snippy_contig_gubbins/
run_gubbins \
    --prefix gubbins \
    --outgroup GCA_000369945.1 \
    --threads 12 \
    ${SUB_PROJECT_DIR}/data/snippy_contig_gubbins/clean.full.aln

docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) -w /data/data/snippy_contig_gubbins snippy:4.6.0_SC2 \
    snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln

## 统计所有序列长度，均为 18366
# awk '/^>/ {if (NR > 1) {print header, length(seq); seq=""} header=$0; next} {seq=seq $0} END {print header, length(seq)}' clean.core.aln
	# 1.	/^>/:
	# •	匹配行首为 > 的标识符行。
	# 2.	if (NR > 1):
	# •	如果当前行不是第一行的 >，输出上一个序列的统计结果。
	# 3.	print header, length(seq):
	# •	输出上一个序列的标识符（> 行内容）和序列长度。
	# 4.	seq="":
	# •	清空 seq 变量，为下一段序列重新累积。
	# 5.	header=$0:
	# •	将当前 > 行的内容存储为 header，用于下一次输出。
	# 6.	seq=seq $0:
	# •	将非 > 的行拼接到变量 seq 中，累积序列数据。
	# 7.	END {print header, length(seq)}:
	# •	处理文件结束时输出最后一个序列的标识符和长度。

set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"