library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

# load('./r_image/compare_lineage_yunnan_v2.RData')
load('./r_image/metadata_v3_with_pairwise.RData')


p_mrm_results   <- readRDS('./r_image/lineage_core_genes_p_mrm_results_panaroo.RDS')
ani_mrm_results <- readRDS('./r_image/lineage_core_genes_ani_mrm_results_panaroo.RDS')

# dnds_results <- readRDS('./r_image/core_genes_codeml_m1a_m2a.RDS')


# protein_distances <- readRDS('./r_image/lineage_core_genes_protein_distances_non_gap.RDS')
# core_genes <- colnames(protein_distances)
# core_genes <- core_genes[6:length(core_genes)]
# rm(protein_distances)

core_genes <- names(ani_mrm_results)

# predictor_cols <- c(
#   "distance_km",
#   "time_interval",
#   "yearly_surface_pressure_2m_difference",
#   "yearly_dewpoint_temperature_2m_difference",
#   "yearly_temperature_2m_difference",
#   "yearly_total_precipitation_sum_difference"
# )

predictor_cols <- c(
  "distance_km"
)

# 构造系数矩阵
mrm_coef_p <- do.call(rbind, lapply(p_mrm_results, function(res) if (!is.null(res)) res$coef else rep(NA, length(predictor_cols))))
rownames(mrm_coef_p) <- core_genes
colnames(mrm_coef_p) <- predictor_cols
# 构造 p 值矩阵
mrm_pvalue_p <- do.call(rbind, lapply(p_mrm_results, function(res) if (!is.null(res)) res$pval else rep(NA, length(predictor_cols))))
rownames(mrm_pvalue_p) <- core_genes
colnames(mrm_pvalue_p) <- predictor_cols

# sum(mrm_pvalue_p[dnds_results$gene, 'distance_km'] < 0.01, na.rm = TRUE)
# 24
# sum(mrm_pvalue_p[setdiff(rownames(mrm_pvalue_p), dnds_results$gene), 'distance_km'] < 0.01, na.rm = TRUE)
# 6

# 原理等同 103 个样本中至少存在 3 个不同的序列
# mrm_coef_p <- as.data.frame(mrm_coef_p)[dnds_results$gene, ]
# mrm_pvalue_p <- as.data.frame(mrm_pvalue_p)[dnds_results$gene, ]

mrm_coef_p <- as.data.frame(mrm_coef_p)
mrm_pvalue_p <- as.data.frame(mrm_pvalue_p)

# 构造 q 值矩阵
mrm_qvalue_p <- apply(mrm_pvalue_p, 2, function(p) p.adjust(p, method = "fdr"))
mrm_qvalue_p <- as.data.frame(mrm_qvalue_p)



# 构造系数矩阵
mrm_coef_ani <- do.call(rbind, lapply(ani_mrm_results, function(res) if (!is.null(res)) res$coef else rep(NA, length(predictor_cols))))
rownames(mrm_coef_ani) <- core_genes
colnames(mrm_coef_ani) <- predictor_cols
# 构造 p 值矩阵
mrm_pvalue_ani <- do.call(rbind, lapply(ani_mrm_results, function(res) if (!is.null(res)) res$pval else rep(NA, length(predictor_cols))))
rownames(mrm_pvalue_ani) <- core_genes
colnames(mrm_pvalue_ani) <- predictor_cols

# 原理等同 103 个样本中至少存在 3 个不同的序列
# mrm_coef_ani <- as.data.frame(mrm_coef_ani)[dnds_results$gene, ]
# mrm_pvalue_ani <- as.data.frame(mrm_pvalue_ani)[dnds_results$gene, ]

mrm_coef_ani <- as.data.frame(mrm_coef_ani)
mrm_pvalue_ani <- as.data.frame(mrm_pvalue_ani)

mrm_pvalue_ani <- rownames_to_column(mrm_pvalue_ani, var = "gene")
mrm_pvalue_ani <- as.data.frame(mrm_pvalue_ani[mrm_pvalue_ani$distance_km < 1 & !is.na(mrm_pvalue_ani$distance_km), ])

mrm_coef_ani <- rownames_to_column(mrm_coef_ani, var = "gene")
mrm_coef_ani <- dplyr::filter(mrm_coef_ani, gene %in% mrm_pvalue_ani$gene)

# 构造 q 值矩阵
mrm_qvalue_ani <- p.adjust(mrm_pvalue_ani$distance_km, method = "fdr")
mrm_qvalue_ani <- as.data.frame(mrm_qvalue_ani)

plot(mrm_coef_p$distance_km, mrm_coef_ani$distance_km)

plot(-log10(as.data.frame(mrm_pvalue_ani)$distance_km), as.data.frame(mrm_coef_ani)$distance_km)



load('./r_image/enricher_db_panaroo_pan.RData')

gene_list <- mrm_coef_p$distance_km
names(gene_list) <- rownames(mrm_coef_p)
gene_list <- na.omit(gene_list)
# gene_list <- gene_list[gene_list > 0]
# gene_list <- abs(gene_list)
# gene_list <- -gene_list
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- clusterProfiler::GSEA(
  gene = gene_list,
  TERM2GENE = db.GO.emapper$TERM2GENE,
  TERM2NAME = db.GO.emapper$TERM2NAME, 
  pvalueCutoff = 1,
  pAdjustMethod = "BH"
)

view(gsea_res@result)




