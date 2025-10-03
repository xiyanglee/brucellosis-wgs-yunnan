library(tidyverse)
library(Biostrings)
library(future.apply)

# 设置工作目录
setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

# 加载样本及核心基因信息
load('./r_image/metadata_v3_with_pairwise.RData')



normalize_datamat <- function(datmat) {
  colnames(datmat) <- datmat[1, ]
  datmat <- datmat[-1, ]
  
  rownames_datamat <- datmat$Gene
  datmat$Gene <- NULL
  
  datmat <- datmat %>% reframe(across(everything(), as.integer))
  rownames(datmat) <- rownames_datamat
  datmat <- as.data.frame(t(datmat))
  
  return(datmat)
}

datmat.onehot <- read.table('./data/panaroo_result_strict_pan/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat.onehot <- normalize_datamat(datmat.onehot)

tip_11  <- metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == '1.1']
tip_13  <- metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == '1.3']
rate_13 <- colMeans(datmat.onehot[tip_13, ])
rate_11 <- colMeans(datmat.onehot[tip_11, ])
# 获取 core genes
core_genes <- intersect(names(rate_11[rate_11 == 1]), names(rate_13[rate_13 == 1]))

rm(tip_11, tip_13, rate_11, rate_13)

# 加载基因存在矩阵（Roary 输出）
# gene_presence_absence <- read.csv('./data/roary_result/gene_presence_absence.csv')
# gene_presence_absence <- gene_presence_absence[, c('Gene', 'Non.unique.Gene.name', 'Annotation', paste0('X', metadata.sample$tip))]
# 
# datmat_name <- gene_presence_absence[, c('Gene', paste0('X', metadata.sample$assembly_accession))]
# datmat_name <- column_to_rownames(datmat_name, var = 'Gene')
# datmat_name <- as.data.frame(t(datmat_name))
# rownames(datmat_name) <- gsub('^X', '', rownames(datmat_name))

# 检查比对文件是否缺失
non_gene_list <- c()
for (GENE_NAME in core_genes) {
  file_path <- paste0('./data/panaroo_result/aligned_gene_sequences/', GENE_NAME, '.aln.fas')
  if (!file.exists(file_path)) {
    message("Missing file: ", file_path)
    non_gene_list <- c(non_gene_list, GENE_NAME)
  }
}; rm(GENE_NAME, file_path)

aligned_genes <- setdiff(core_genes, non_gene_list)

# 并行计算每个基因的 p-distance 和 1 - ANI
plan(multisession)

dist_matrix_list <- future_lapply(aligned_genes, function(GENE_NAME) {
  file_path <- paste0('./data/panaroo_result/aligned_gene_sequences/', GENE_NAME, '.aln.fas')
  if (!file.exists(file_path)) return(NULL)
  
  fasta_file <- readDNAStringSet(file_path)
  
  # 删除分号后的部分
  names(fasta_file) <- sub(";.*", "", names(fasta_file))
  
  # name_map <- setNames(rownames(datmat_name), datmat_name[, GENE_NAME])
  # matched_names <- intersect(names(fasta_file), names(name_map))
  # if (length(matched_names) < 2) return(NULL)
  
  matched_names <- intersect(names(fasta_file), metadata.sample$assembly_accession)
  aligned_seqs <- fasta_file[matched_names]
  # names(aligned_seqs)[names(aligned_seqs) %in% matched_names] <- name_map[matched_names]
  
  seq_matrix <- as.matrix(aligned_seqs)
  seq_names <- names(aligned_seqs)
  
  n <- length(seq_names)
  p_mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(seq_names, seq_names))
  ani_mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(seq_names, seq_names))
  
  for (i in seq_along(seq_names)) {
    for (j in seq_along(seq_names)) {
      seq1 <- seq_matrix[i, ]
      seq2 <- seq_matrix[j, ]
      
      valid_pos <- which(seq1 != "-" & seq2 != "-")
      p_dist <- if (length(valid_pos) == 0) NA else {
        sum(seq1[valid_pos] != seq2[valid_pos]) / length(valid_pos)
      }
      
      ani_sim <- sum(seq1 == seq2) / length(seq1)
      ani_diff <- 1 - ani_sim
      
      p_mat[i, j] <- p_dist
      ani_mat[i, j] <- ani_diff
    }
  }
  
  return(list(p = p_mat, ani = ani_mat))
})

names(dist_matrix_list) <- aligned_genes

# 整合样本对结构
# gene_list <- names(dist_matrix_list)
# codon_distances <- pairwise_metadata.sample %>%
#   dplyr::select(tip_pair, tip_1, tip_2, lineage_1, lineage_2)


# 保存输出文件
# saveRDS(dist_matrix_list, './r_image/lineage_core_genes_codon_distances_p_and_ani_panaroo.RDS')
