library(fastbaps)
library(ape)
library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/clustering_analysis_fastBAPS.RData')



# 合并分类群
partition.merge <- data.frame(tip = raxml_tree$tip.label) %>% 
  left_join(partition.baps, by = "tip") %>% 
  left_join(partition.opt, by = "tip") %>% 
  left_join(partition.tree, by = "tip")


## 查看 cluster 及其具体位置
## 便于确定分类群位置
# ggtree(raxml_tree) %<+% partition.tree + geom_tippoint(aes(color = group.tree), size = 2) + 
#   paletteer::scale_color_paletteer_d("ggsci::category20_d3")
# 
# ggtree(raxml_tree) %<+% 
#   (
#     partition.merge %>% 
#       mutate(marked = case_when(
#         group.tree == 19 ~ "label",
#         .default = NA
#       ))
#   ) + 
#   geom_tippoint(aes(color = marked), size = 2) + 
#   ggsci::scale_color_bmj()

partition.merge$lineage_level_2a <- NA
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(11, 12, 13, 14, 15)] <- '1.1'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(5, 6, 7, 8, 9)] <- '1.2'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(16)] <- '1.3a'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(10)] <- '1.3b'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(19)] <- '1.u'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(17)] <- '2.1'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(3)] <- '2.2'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(18)] <- '2.3'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(1)] <- '3.1'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(4)] <- '3.2'
partition.merge$lineage_level_2a[partition.merge$group.tree %in% c(2)] <- '3.3'

partition.merge <- partition.merge %>% 
  mutate(
    lineage_level_2 = case_match(
      lineage_level_2a, 
      c('1.3a', '1.3b') ~ '1.3', 
      .default = lineage_level_2a
    ), 
    lineage_level_1 = case_match(
      lineage_level_2, 
      c('1.1', '1.2', '1.3', '1.u') ~ '1', 
      c('2.1', '2.2', '2.3') ~ '2', 
      c('3.1', '3.2', '3.3') ~ '3', 
      .default = NA
    )
  )

p <- ggtree(raxml_tree)
p <- facet_plot(p, panel = "fastBAPS unoptimized\nInitial clusters based on SNPs", data = partition.merge, geom = geom_tile, aes(x = as.integer(group.baps)), color = '#7f9faf')
p <- facet_plot(p, panel = "fastBAPS optimized\nInitial clusters based on SNPs", data = partition.merge, geom = geom_tile, aes(x = as.integer(group.opt)), color = '#cf929e')
p <- facet_plot(p, panel = "fastBAPS optimized\nInitial clusters based on RAxML tree", data = partition.merge, geom = geom_tile, aes(x = as.integer(group.tree)), color = '#06436f')
p <- facet_plot(p, panel = "ANI", data = partition.merge, geom = geom_tile, aes(x = as.integer(as.factor(lineage_level_2))), color = '#f8b862')
p

# ggsave(filename = paste0("./plot/clustering_analysis/clusters_by_diff_methods.pdf"), height = 6, width = 8)

saveRDS(partition.merge, './r_image/clustering_analysis_designating_lineages.RDS')
