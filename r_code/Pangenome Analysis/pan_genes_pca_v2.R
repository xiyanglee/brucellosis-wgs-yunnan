# install.packages("ggalt", lib = "/public/r_share_library/4.1")
library(tidyverse)
library(ggalt)
library(plotly)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/core_genes.RData')
rm(resample_res, subsample)
load('./r_image/compare_lineage_yunnan_v2.RData')



datmat <- read.table('./data/roary_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat <- t(datmat) %>% as.data.frame()

colnames(datmat) <- datmat[1, ]
datmat <- datmat[-1, ]

rownames_datamat <- datmat$Gene
datmat$Gene <- NULL
datmat <- datmat %>% reframe(across(everything(), as.integer))
rownames(datmat) <- rownames_datamat
rm(rownames_datamat)

pca_datmat <- datmat %>%
  select(where(~ {
    mean_val <- mean(.)
    mean_val >= 0.002 # mean_val <= 0.9 &
  })) %>%
  select(where(~ sd(.) > 0)) %>% 
  prcomp(., scale. = FALSE)

pca_data <- data.frame(
  tip = rownames(pca_datmat$x),
  PC1 = pca_datmat$x[, 1], 
  PC2 = pca_datmat$x[, 2] 
) %>% 
  right_join(., metadata_tree, by = "tip") %>% 
  dplyr::filter(tip != "Reference")

pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = group_2)) + # continent area geographic_country
  geom_point(size = 1.5, alpha = 0.5) + 
  geom_encircle(data = pca_data %>% filter(group_2 == "1.3"), # 筛选目标类别
                aes(x = PC1, y = PC2),
                linetype = 'dashed',
                color = "#7C5CA2", 
                size = 2, 
                s_shape = 2,
                expand = 0.05) + # 设置边距
  # labs(
  #   title = "PCA of dotmat with Metadata",
  #   x = "Principal Component 1",
  #   y = "Principal Component 2",
  #   color = "Group"
  # ) +
  theme_void(base_size = 14) + 
  # theme(
  #   plot.margin = margin(10, 10, 10, 10) # 调整边距
  # ) +
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

# ggsave('./plot/pan_genome_analysis/pan_genes_pca_col2.pdf', width = 4, height = 3)

################################################################################
## PCoA 分析

library(vegan) # 用于 PCoA 分析

# 1. 计算距离矩阵
dist_matrix <- datmat %>% 
  # select(where(~ {
  #   mean_val <- mean(.)
  #   mean_val >= 0.001
  # })) %>%
  vegdist(., method = "euclidean") # euclidean bray

# 2. 执行 PCoA 分析
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)  # k = 2 表示提取 2 个维度

# 3. 提取 PCoA 结果并合并到 metadata
pcoa_scores <- data.frame(
  tip = rownames(pcoa_result$points),
  PC1 = pcoa_result$points[, 1], 
  PC2 = pcoa_result$points[, 2] 
) %>% 
  right_join(., metadata_tree, by = "tip") %>% 
  dplyr::filter(tip != "Reference")

# {
#   ggplot(pcoa_scores, aes(x = PC1, y = PC2, color = group_2)) + # continent collection_year area geographic_country
#     geom_point(size = 2) +
#     geom_text(nudge_y = 0.1, size = 3) +  # 添加样本标签
#     theme_minimal()
#   } %>% 
#   ggplotly(tooltip = "color")

pcoa_scores %>% 
  mutate(PC2 = -PC2) %>% 
  ggplot(aes(x = PC1, y = PC2, color = group_2)) + # continent area geographic_country
  geom_point(size = 1.5, color = "black", fill = NA, alpha = 1, shape = 21, stroke = 1.2) + 
  geom_point(size = 1.5, alpha = 1) + 
  geom_point(data = dplyr::filter(pcoa_scores, group_2 == "1.3") %>% mutate(PC2 = -PC2), size = 1.5, alpha = 1) + 
  geom_encircle(data = pca_data %>% filter(group_2 == "1.3"), # 筛选目标类别
                aes(x = PC1, y = PC2),
                linetype = 'dashed',
                color = "#632fa1", 
                size = 3, 
                s_shape = 1.2,
                expand = 0.05) + # 设置边距
  theme_void(base_size = 14) + 
  annotate("segment", 
           x = min(pca_data$PC1) - 1, xend = min(pca_data$PC1) + 2, 
           y = min(pca_data$PC2) - 1, yend = min(pca_data$PC2) - 1, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + # x 轴箭头
  annotate("segment", 
           x = min(pca_data$PC1) - 1, xend = min(pca_data$PC1) - 1, 
           y = min(pca_data$PC2) - 1, yend = min(pca_data$PC2) + 2, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + # y 轴箭头
  annotate("text", x = min(pca_data$PC1) + 2.5, y = min(pca_data$PC2) - 2.5, 
           label = "PCo1", hjust = 1, vjust = 1, color = "black") + # x 轴标签
  annotate("text", x = min(pca_data$PC1) - 2.5, y = min(pca_data$PC2) + 2.5, 
           label = "PCo2", hjust = 1, vjust = 0, color = "black", angle = 90) + # y 轴标签 
  coord_equal() + 
  scale_color_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                                "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                                "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                     na.value = "transparent")

# ggsave('./plot/pan_genome_analysis/pan_genes_pcoa_col2_type_2.pdf', width = 4, height = 4)
