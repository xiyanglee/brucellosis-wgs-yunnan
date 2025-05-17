#!/bin/bash

## chmod +x 12_2_all_genes_codeml_m1a_m2a_m7_m8.sh
## ./12_2_all_genes_codeml_m1a_m2a_m7_m8.sh

# 设置最大并发数
MAX_JOBS=10

# 根目录
BASE_DIR="/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly/data/dnds_per_gene_m1a_vs_m2a_m7_m8_all_genes"

# 导出执行函数以供 parallel 使用
run_codeml() {
    local dir="$1"
    cd "$dir" || exit 1
    echo "Running codeml in $(basename "$dir")"
    
    codeml codeml.m1a.ctl > m1a.log 2>&1
    codeml codeml.m2a.ctl > m2a.log 2>&1
    codeml codeml.m7.ctl > m7.log 2>&1
    codeml codeml.m8.ctl > m8.log 2>&1
}
export -f run_codeml

# 使用 find 获取所有子目录（假设每个目录都有 codeml 文件）
find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d \
    | parallel -j "$MAX_JOBS" run_codeml {}

