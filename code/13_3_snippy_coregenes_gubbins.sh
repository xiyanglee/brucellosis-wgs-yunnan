#!/bin/bash

printf "############# JOB INFO START #############\n"
printf "%15s %s\n" "START_TIME|"     "`date`"
printf "##########################################\n"

set -uex

source _data_dir.sh

# nohup bash 13_3_snippy_coregenes_gubbins.sh >> _log/snippy_coregenes_gubbins.log 2>&1 &


docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) snippy:4.6.0_SC2 \
bash -c "snippy-clean_full_aln /data/data/coregene_merge_snippy_alignment_merge/core.full.aln > /data/data/coregene_merge_snippy_gubbins/clean.full.aln"

# --filter_percentage FILTER_PERCENTAGE, -f FILTER_PERCENTAGE Filter out taxa with more than this percentage of gaps (default: 25)
# 默认参数下以下序列会被过滤掉，包括 outgroup
# Filtering input alignment...
# Excluded sequence GCA_000369945.1 because it had 27.91424525843529 percentage missing data while a maximum of 25 is allowed
# Excluded sequence GCA_947241995.1 because it had 39.34085892883493 percentage missing data while a maximum of 25 is allowed
# Excluded sequence GCA_947242235.1 because it had 43.46219777528704 percentage missing data while a maximum of 25 is allowed
# Excluded sequence GCA_947242715.1 because it had 59.72927321647782 percentage missing data while a maximum of 25 is allowed
cd ${SUB_PROJECT_DIR}/data/coregene_merge_snippy_gubbins/
  run_gubbins \
--prefix gubbins \
--outgroup GCA_000369945.1 \
--filter_percentage 60 \
--threads 12 \
${SUB_PROJECT_DIR}/data/coregene_merge_snippy_gubbins/clean.full.aln

docker run --rm -v ${SUB_PROJECT_DIR}:/data -u $(id -u):$(id -g) -w /data/data/coregene_merge_snippy_gubbins snippy:4.6.0_SC2 \
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln

set +x
printf "############## JOB INFO END ##############\n"
printf "%15s %s\n" "END_TIME|"     "`date`"
printf "##########################################\n"