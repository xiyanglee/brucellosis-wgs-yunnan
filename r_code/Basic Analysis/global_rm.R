library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/compare_lineage_yunnan_v2.RData')


gubbins_results <- read.table('./data/snippy_contig_gubbins/gubbins.per_branch_statistics.csv', sep = '\t', header = TRUE)
gubbins_results <- inner_join(gubbins_results, metadata_tree, by = c("Node" = "tip"))

table(gubbins_results$r.m)

gubbins_summary <- gubbins_results %>% 
  dplyr::select(r.m, group_2) %>% 
  group_by(group_2) %>% 
  summarise(r.m = mean(r.m))

# ggplot() + 
#   geom_jitter(data = gubbins_results, aes(x = group_2, y = r.m, fill = group_2), color = "black", shape = 21, size = 2) + 
#   geom_col(data = gubbins_summary, aes(x = group_2, y = r.m, fill = group_2), alpha = 0.5, color = 'black') + 
#   scale_y_sqrt(expand = c(0, 0, 0, 0.1)) + 
#   theme_bw(base_size = 14)

ggplot() + 
  geom_jitter(data = gubbins_results, aes(x = group_2, y = r.m, fill = group_2), 
              color = "black", shape = 21, size = 2) + 
  geom_col(data = gubbins_summary, aes(x = group_2, y = r.m, fill = group_2), 
           alpha = 0.75, color = 'black') + 
  scale_y_sqrt(expand = c(0, 0, 0, 0.1)) + 
  labs(x = "Lineage", y = "r/m") +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_blank(),            # 去除整个 panel border
    axis.line = element_line(color = "black"), # 添加坐标轴线（默认只会显示左和下）
    legend.position = "none",                  # 去除图例
    panel.grid.minor.x = element_blank(),      # 去除垂直方向的次要 grid
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"), 
    panel.grid.minor.y = element_blank()
  ) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                               "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                               "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                    na.value = "transparent")


# ggsave('./plot/pan_genome_analysis/lineage_rm_global_col2.pdf', width = 4, height = 4)
