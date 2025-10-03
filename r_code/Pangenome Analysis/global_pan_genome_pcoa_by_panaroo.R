library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



datmat <- read.table('./data/panaroo_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat <- t(datmat) %>% as.data.frame()

colnames(datmat) <- datmat[1, ]
datmat <- datmat[-1, ]

rownames_datamat <- datmat$Gene
datmat$Gene <- NULL
datmat <- datmat %>% reframe(across(everything(), as.integer))
rownames(datmat) <- rownames_datamat
rm(rownames_datamat)

metadata.phylo <- metadata.phylo %>% dplyr::filter(lineage_level_2 != "1.u")
datmat <- datmat[metadata.phylo$assembly_accession, , drop = FALSE]



################################################################################
## PCA 分析

pca_datmat <- prcomp(datmat, scale. = FALSE)
pca_data   <- data.frame(
  assembly_accession = rownames(pca_datmat$x),
  PC1 = pca_datmat$x[, 1], 
  PC2 = pca_datmat$x[, 2] 
) %>% 
  right_join(., metadata.phylo, by = "assembly_accession")

ggplot(pca_data, aes(x = PC1, y = PC2, color = lineage_level_2)) + # continent area geographic_country
  geom_point(size = 1.5, alpha = 0.5) + 
  # labs(
  #   title = "PCA of dotmat with Metadata",
  #   x = "Principal Component 1",
  #   y = "Principal Component 2",
  #   color = "Group"
  # ) +
  theme_void(base_size = 14) +
  theme(
    plot.margin = margin(10, 10, 10, 10) # 调整边距
  ) +
  # 添加左下角的箭头
  annotate("segment", 
           x = min(pca_data$PC1) - 1, xend = min(pca_data$PC1) + 2, 
           y = min(pca_data$PC2) - 1, yend = min(pca_data$PC2) - 1, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + # x 轴箭头
  annotate("segment", 
           x = min(pca_data$PC1) - 1, xend = min(pca_data$PC1) - 1, 
           y = min(pca_data$PC2) - 1, yend = min(pca_data$PC2) + 2, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + # y 轴箭头
  annotate("text", x = min(pca_data$PC1) + 2.5, y = min(pca_data$PC2) - 2.5, 
           label = "PC1", hjust = 1, vjust = 1, color = "black") + # x 轴标签
  annotate("text", x = min(pca_data$PC1) - 2.5, y = min(pca_data$PC2) + 2.5, 
           label = "PC2", hjust = 1, vjust = 0, color = "black", angle = 90) + # y 轴标签 
  coord_equal() + 
  scale_color_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E",
                                "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5",
                                "3.1" = "#c9a7b7", "3.2" = "#EBCACB"),
                     na.value = "transparent")

# ggsave('./plot/pan_genome_analysis/global_pan_genome_pca_by_panaroo.pdf', width = 4, height = 3)



################################################################################
## PCoA 分析

library(vegan) # 用于 PCoA 分析

# 1. 计算距离矩阵
dist_matrix <- vegdist(datmat, method = "kulczynski") # manhattan jaccard kulczynski
# 2. 执行 PCoA 分析
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)  # k = 2 表示提取 2 个维度

# 3. 提取 PCoA 结果并合并到 metadata
pcoa_scores <- data.frame(
  assembly_accession = rownames(pcoa_result$points),
  PC1 = pcoa_result$points[, 1], 
  PC2 = pcoa_result$points[, 2] 
) %>% 
  right_join(., metadata.phylo, by = "assembly_accession")

pcoa_scores %>% 
  ggplot(aes(x = PC1, y = PC2, color = lineage_level_2)) + # continent area geographic_country
  geom_point(size = 2, color = "black", fill = NA, alpha = 1, shape = 21, stroke = 1.2) + 
  geom_point(size = 2, alpha = 1) + 
  geom_point(data = dplyr::filter(pcoa_scores, lineage_level_1 == "3"), size = 2, alpha = 1) + 
  geom_point(data = dplyr::filter(pcoa_scores, lineage_level_2 == "1.3"), size = 2, alpha = 1) + 
  stat_ellipse(data = pcoa_scores %>% filter(lineage_level_1 == "1"),
               aes(x = PC1, y = PC2),
               level = 0.95,
               linetype = 'dashed',
               color = "#2c9969",
               linewidth = 0.75,
               type = "norm") + 
  stat_ellipse(data = pcoa_scores %>% filter(lineage_level_1 == "2"),
               aes(x = PC1, y = PC2),
               level = 0.95,
               linetype = 'dashed',
               color = "#c4901d",
               linewidth = 0.75,
               type = "norm") + 
  stat_ellipse(data = pcoa_scores %>% filter(lineage_level_1 == "3"),
               aes(x = PC1, y = PC2),
               level = 0.95,
               linetype = 'dashed',
               color = "#b86c90",
               linewidth = 0.75,
               type = "norm") + 
  theme_void(base_size = 14) + 
  annotate("segment",
           x = min(pcoa_scores$PC1) - 0.01, xend = min(pcoa_scores$PC1) + 0.01,
           y = min(pcoa_scores$PC2) - 0.01, yend = min(pcoa_scores$PC2) - 0.01,
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + # x 轴箭头
  annotate("segment",
           x = min(pcoa_scores$PC1) - 0.01, xend = min(pcoa_scores$PC1) - 0.01,
           y = min(pcoa_scores$PC2) - 0.01, yend = min(pcoa_scores$PC2) + 0.01,
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + # y 轴箭头
  annotate("text", x = min(pcoa_scores$PC1) + 0.01, y = min(pcoa_scores$PC2) - 0.01,
           label = "PCo1", hjust = 0.01, vjust = 0.01, color = "black") + # x 轴标签
  annotate("text", x = min(pcoa_scores$PC1) - 0.01, y = min(pcoa_scores$PC2) + 0.01,
           label = "PCo2", hjust = 0.01, vjust = 0, color = "black", angle = 90) + # y 轴标签
  coord_equal() +
  
  scale_color_manual(name = "Lineage", 
                     values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                                "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                                "3.1" = "#c9a7b7", "3.2" = "#EBCACB", "3.3" = "#EBCACB"), 
                     na.value = "black")

# ggsave('./plot/pan_genome_analysis/global_pan_genome_pcoa_by_panaroo.pdf', width = 4, height = 4)



################################################################################
## PERMANOVA 分析

dist_matrix   <- vegdist(datmat, method = "jaccard")
adonis_level_1 <- adonis2(dist_matrix ~ lineage_level_1, data = metadata.phylo, permutations = 999)

print(adonis_level_1)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_matrix ~ lineage_level_1, data = metadata.phylo, permutations = 999)
#                  Df SumOfSqs      R2      F Pr(>F)    
# lineage_level_1   2  0.19335 0.10751 50.776  0.001 ***
# Residual        843  1.60502 0.89249                  
# Total           845  1.79836 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis_level_2 <- adonis2(dist_matrix ~ lineage_level_2, data = metadata.phylo, permutations = 999)
print(adonis_level_2)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_matrix ~ lineage_level_2, data = metadata.phylo, permutations = 999)
#                  Df SumOfSqs      R2      F Pr(>F)    
# lineage_level_2   8  0.61009 0.33925 53.717  0.001 ***
# Residual        837  1.18827 0.66075                  
# Total           845  1.79836 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



################################################################################
## pairwise PERMANOVA 分析

# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"))
library(pairwiseAdonis)

pairwise_level_1 <- pairwise.adonis2(dist_matrix ~ lineage_level_1, data = metadata.phylo, perm = 999)

# 提取所有组对名（排除非结果列表元素，如 parent_call）
pairwise_keys <- names(pairwise_level_1)[!names(pairwise_level_1) %in% "parent_call"]
# 提取每个组对的统计信息
pairwise_level_1_results <- do.call(rbind, lapply(pairwise_keys, function(key) {
  res <- pairwise_level_1[[key]]
  if (!is.null(res)) {
    stats_row <- res[1, ]  # 仅提取第一行（lineage_level_2 的检验结果）
    data.frame(
      group = key,
      Df = stats_row$Df,
      # SumOfSqs = stats_row$SumOfSqs,
      R2 = stats_row$R2,
      F = stats_row$F,
      p_value = stats_row$`Pr(>F)`,
      row.names = NULL
    )
  }
}))

pairwise_level_2 <- pairwise.adonis2(dist_matrix ~ lineage_level_2, data = metadata.phylo, perm = 999)

pairwise_keys <- names(pairwise_level_2)[!names(pairwise_level_2) %in% "parent_call"]
pairwise_level_2_results <- do.call(rbind, lapply(pairwise_keys, function(key) {
  res <- pairwise_level_2[[key]]
  if (!is.null(res)) {
    stats_row <- res[1, ]  # 仅提取第一行（lineage_level_2 的检验结果）
    data.frame(
      group = key,
      Df = stats_row$Df,
      R2 = stats_row$R2,
      F = stats_row$F,
      p_value = stats_row$`Pr(>F)`,
      row.names = NULL
    )
  }
}))

# write.csv(pairwise_level_1_results, './result/pan_genome_analysis/lineage_level_1_pangenome_pairwise_PERMANOVA_by_panaroo.csv', row.names = FALSE)
# write.csv(pairwise_level_2_results, './result/pan_genome_analysis/lineage_level_2_pangenome_pairwise_PERMANOVA_by_panaroo.csv', row.names = FALSE)

