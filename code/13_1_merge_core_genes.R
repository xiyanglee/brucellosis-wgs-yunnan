library(Biostrings)
library(tidyverse)
# library(ape)

## 以下读取合并的是 core genes (>99%) 的 codon alignment，合并后建树，但 SNP 过多 (120w+)，可能是由于 codon 对齐导致错位较多

# datmat <- read.table('./data/roary_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
# datmat <- t(datmat) %>% as.data.frame()
# colnames(datmat) <- datmat[1, ]
# datmat <- datmat[-1, ]
# # selected_indices <- grepl("^20", datmat$Gene)  # 筛选样本名以 "20" 开头的行
# # datmat <- datmat[selected_indices, ]
# rownames_datamat <- datmat$Gene
# datmat$Gene <- NULL
# datmat <- datmat %>% reframe(across(everything(), as.integer))
# rownames(datmat) <- rownames_datamat
# rm(rownames_datamat)
# 
# coregene_names <- colnames(datmat)[colMeans(datmat) > 0.99]
# sample_names <- rownames(datmat)
# 
# 
# parse_id <- function(id) {
#   # 1) 去掉 ">" 和 "PROKKATAG-"
#   id2 <- sub("^>", "", id)
#   id2 <- sub("^PROKKATAG-", "", id2)
#   
#   # 2) 区分是否以 "GCA_" 开头
#   if (grepl("^GCA_", id2)) {
#     # 2.1) 若以 GCA_ 开头 => 先截取到第二个 '_' 之前
#     #      例如 "GCA_000022625.1_02550" => "GCA_000022625.1"
#     part <- sub("^([^_]+_[^_]+)_.*$", "\\1", id2)
#   } else {
#     # 2.2) 若不以 GCA_ 开头 => 只截取到第一个 '_' 之前
#     #      例如 "201902_03111" => "201902"
#     part <- sub("^([^_]+)_.*$", "\\1", id2)
#   }
#   
#   part
# }
# 
# # 初始化存储每个基因比对结果的列表
# aligned_list <- list()
# 
# for (gene in coregene_names) {
#   print(gene)
#   
#   # 读取每个基因的 codon 比对文件
#   aln <- readDNAStringSet(paste0('./data/pan_cds_align/', gene, '_codon_aligned.fasta'))
#   names(aln) <- sapply(names(aln), parse_id)
#   
#   # 记录比对长度（每条序列长度应该一致）
#   aln_width <- width(aln)[1]
#   
#   # 补齐缺失样本
#   missing_samples <- setdiff(sample_names, names(aln))
#   if (length(missing_samples) > 0) {
#     gap_seqs <- DNAStringSet(rep(paste0(rep("-", aln_width), collapse = ""), length(missing_samples)))
#     names(gap_seqs) <- missing_samples
#     aln <- c(aln, gap_seqs)
#   }
#   
#   # 按 sample_names 顺序排列
#   aln <- aln[sample_names]
#   
#   aligned_list[[gene]] <- aln
# }; rm(gene, aln, aln_width, missing_samples, gap_seqs)
# 
# # 拼接所有基因的比对结果
# # 使用 Reduce + xscat 拼接所有 DNAStringSet
# concatenated_alignment <- Reduce(xscat, aligned_list)
# 
# # 获取样本数和序列长度
# n_samples <- length(concatenated_alignment)
# seq_length <- width(concatenated_alignment)[1]  # 所有序列应等长
# 
# # 创建空向量保存是否为SNP
# is_snp <- logical(seq_length)
# 
# # 将整个对象转换为字符向量（每条是完整序列）
# seqs_char <- as.character(concatenated_alignment)
# 
# # 按列检查多态性（变异在至少两个样本中出现）
# for (i in seq_len(seq_length)) {
#   col_i <- substr(seqs_char, i, i)
#   # 去除不确定碱基
#   col_i <- col_i[col_i != "n" & col_i != "N"]
#   # 统计每种碱基的频次（包括"-"）
#   base_counts <- table(col_i)
#   # 有效碱基类型中，至少两种且至少一种频数 >= 2（即变异存在于多个样本）
#   is_snp[i] <- length(base_counts) > 1 && any(base_counts >= 2)
# }; rm(i, col_i, base_counts)
# 
# # 从原始字符串中提取 SNP 位点，重新组合为字符序列
# snp_seqs <- vapply(seqs_char, function(seq) {
#   paste0(strsplit(seq, "")[[1]][is_snp], collapse = "")
# }, character(1))
# 
# # 转为 DNAbin 格式
# names(snp_seqs) <- names(concatenated_alignment)
# snp_dna <- as.DNAbin(strsplit(snp_seqs, split = ""))
# 
# rm(n_samples, seq_length, is_snp, seqs_char, snp_seqs)
# gc()
# 
# # save.image('./r_image/core_codon_alignment.RData')
# 
# tree <- nj(dist.dna(snp_dna, model = "N", pairwise.deletion = TRUE))




# 读取 gene_presence_absence.Rtab 表格，提取 core genes 和样本名
datmat <- read.table('./data/roary_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat <- t(datmat) %>% as.data.frame()
colnames(datmat) <- datmat[1, ]
datmat <- datmat[-1, ]
rownames_datamat <- datmat$Gene
datmat$Gene <- NULL
datmat <- datmat %>% reframe(across(everything(), as.integer))
rownames(datmat) <- rownames_datamat
rm(rownames_datamat)

coregene_names <- colnames(datmat)[colMeans(datmat) > 0.99]
sample_names <- rownames(datmat)

# 解析 ID 成为样本名
parse_id <- function(id) {
  id2 <- sub("^>", "", id)
  id2 <- sub("^PROKKATAG-", "", id2)
  if (grepl("^GCA_", id2)) {
    part <- sub("^([^_]+_[^_]+)_.*$", "\\1", id2)
  } else {
    part <- sub("^([^_]+)_.*$", "\\1", id2)
  }
  part
}

# 初始化样本合并序列列表
sample_concat_seqs <- setNames(vector("list", length(sample_names)), sample_names)
for (s in sample_names) sample_concat_seqs[[s]] <- DNAString()

# 遍历每个核心基因，读取 _cds.fasta 文件（非对齐）
for (gene in coregene_names) {
  cat("Processing:", gene, "\n")
  cds_path <- paste0('./data/pan_cds_align/', gene, '_cds.fasta')
  if (!file.exists(cds_path)) next
  
  cds <- readDNAStringSet(cds_path)
  names(cds) <- sapply(names(cds), parse_id)
  
  # 遍历每个样本，拼接其序列（若该样本没有该基因则跳过）
  for (s in intersect(sample_names, names(cds))) {
    sample_concat_seqs[[s]] <- xscat(sample_concat_seqs[[s]], cds[[s]])
  }
}; rm(cds, cds_path, gene)

# 输出每个样本的拼接序列到 fasta 文件
output_dir <- './data/core_genes_fasta/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 保存序列
# for (s in sample_names) {
#   seq <- sample_concat_seqs[[s]]
#   dss <- DNAStringSet(seq)
#   names(dss) <- s
#   writeXStringSet(dss, filepath = file.path(output_dir, paste0(s, '.fasta')))
# }; rm(seq, dss, output_dir, s)

# 获取所有拼接序列长度
seq_lengths <- sapply(sample_concat_seqs, nchar)

# 构建 data.frame
length_df <- data.frame(
  sample = names(seq_lengths),
  length = as.integer(seq_lengths)
)

# saveRDS(length_df, './r_image/core_gene_sequences_merge_length_in_all_samples.RDS')


