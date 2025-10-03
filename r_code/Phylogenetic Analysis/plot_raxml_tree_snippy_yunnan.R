library(tidyverse)
library(ggtree)
library(ape)

load('./r_image/cluster_tree.RData')
load('./r_image/metadata.RData')

metadata_tree <- left_join(group_data_2, metadata, by = c('tip' = 'assembly_accession'))
# metadata_tree <- metadata_tree %>% 
#   mutate(
#     label_yunnan = case_when(
#       continent == 'Yunnan' ~ 'Yunnan', 
#       .default = NA), 
#     continent = case_when(
#       continent == 'Yunnan' ~ 'Asia', 
#       .default = continent
#     )
#   )

metadata_tree.yunnan <- dplyr::filter(metadata_tree, geographic_country == 'Yunnan')
drop_tips <- setdiff(raxml_tree$tip.label, metadata_tree.yunnan$tip)
raxml_tree.yunnan <- drop.tip(raxml_tree, drop_tips)


## 寻找各分类最近祖先节点
phytools::findMRCA(raxml_tree.yunnan, tip = metadata_tree.yunnan$tip[metadata_tree.yunnan$group_2 == '1.1']) # 113
phytools::findMRCA(raxml_tree.yunnan, tip = metadata_tree.yunnan$tip[metadata_tree.yunnan$group_2 == '1.2']) # 191
phytools::findMRCA(raxml_tree.yunnan, tip = metadata_tree.yunnan$tip[metadata_tree.yunnan$group_2 == '1.3']) # 104

ggtree(raxml_tree.yunnan) %<+% metadata_tree.yunnan + 
  geom_hilight(node = 113, alpha = 0.3, fill = '#FBCBD2') +
  geom_hilight(node = 191, alpha = 0.3, fill = '#ee827c') +
  geom_hilight(node = 104, alpha = 0.3, fill = '#c9171e') +
  geom_cladelabel(node = 113, label = "1.1", align = T, color = '#FBCBD2', linewidth = 2) +
  geom_cladelabel(node = 191, label = "1.2", align = T, color = '#ee827c', linewidth = 2) + 
  geom_cladelabel(node = 104, label = "1.3", align = T, color = '#c9171e', linewidth = 2)
   
# ggsave('./plot/tree/raxml_tree_snippy_yunnan.pdf', width = 4, height = 4)



