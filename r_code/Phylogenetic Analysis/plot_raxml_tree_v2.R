# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("gtonkinhill/fastbaps", build_vignettes = TRUE))
# withr::with_libpaths("/public/r_share_library/4.1", remotes::install_github('YuLab-SMU/ggtree'))
# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("xiangpin/ggtreeExtra"))
# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("GuangchuangYu/ggimage"))
# sudo apt-get install libmagick++-dev
# install.packages("ggimage", lib = "/public/r_share_library/4.1")
library(fastbaps)
library(ape)
library(Matrix)
library(tidyverse)
library(ggtree)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

## based on: plot_raxml_tree_snippy.R
load('./r_image/cluster_tree.RData')
load('./r_image/metadata.RData')

metadata_tree <- left_join(group_data_2, metadata, by = c('tip' = 'assembly_accession'))
metadata_tree <- metadata_tree %>% 
  mutate(
    label_yunnan = case_when(
      continent == 'Yunnan' ~ 'Yunnan', 
      .default = NA), 
    continent = case_when(
      continent == 'Yunnan' ~ 'Asia', 
      .default = continent
    )
  )

library(ggtreeExtra)
# ggtree(ape::compute.brlen(raxml_tree, method = "Grafen"), size = 0.3, layout = "fan", open.angle = 90) %<+%

# ggtree(raxml_tree, layout = "fan", open.angle = 90, linewidth = 0.25) %<+%
#   (
#     metadata_tree %>% 
#       mutate(
#         plot_year = case_when(
#           collection_year < 1990 ~ 'before 1990s', 
#           collection_year >= 1990 & collection_year < 2000 ~ '1990s', 
#           collection_year >= 2000 & collection_year < 2010 ~ '2000s', 
#           collection_year >= 2010 & collection_year < 2020 ~ '2010s', 
#           collection_year >= 2020 ~ '2020s', 
#           .default = NA
#         ), 
#         area = case_when(
#           area == "Yunnan" ~ "East Asia", 
#           .default = area
#         )
#       )
#   ) + 
#   geom_tiplab(aes(color = label_yunnan), as_ylab = FALSE, align = T, linesize = 0.5, linetype = 1, offset = 0.001) + 
#   scale_color_manual(values = c("Yunnan" = "#dcdddd"), na.value = "transparent") + 
#   theme(
#     axis.text.y = element_blank()  # 移除 Y 轴刻度文本
#   ) +
#   ggnewscale::new_scale_color() + 
#   geom_tippoint(aes(color = group_2), size = 0.5) + 
#   scale_color_manual(values = c("1.1" = "#FBCBD2", "1.2" = "#ee827c", "1.3" = "#c9171e", 
#                                 "2.1" = "#aed0ee", "2.2" = "#2ca9e1", "2.3" = "#2e59a7", 
#                                 "3.1" = "#c0d695", "3.2" = "#4f6f46"), 
#                      na.value = "transparent") + 
#   ggnewscale::new_scale_color() + 
#   geom_fruit(
#     geom = geom_tile,
#     mapping = aes(y = tip, fill = continent),
#     width = 0.0025,
#     offset = 0.05
#   ) + 
#   scale_fill_manual(values = c("Asia" = "#8ECFC9", "Europe" = "#FFBE7A", "Africa" = "#bba1cb", "America" = "#FA7F6F"),
#                     na.value = "grey95") +
#   ggnewscale::new_scale_fill() + 
#   geom_fruit(
#     geom = geom_tile,
#     mapping = aes(y = tip, fill = area),
#     width = 0.0025,
#     offset = 0.06
#   ) + 
#   scale_fill_manual(values = c("East Asia" = "#4f6f46", 
#                                "Western Asia" = "#99bcac", 
#                                "Southeast Asia" = "#a9be7b", 
#                                "South Asia" = "#c0d695", 
#                                "North Asia" = "#b7d332", 
#                                "Eastern Europe" = "#fce2c4", 
#                                "Central Europe" = "#fedc5e", 
#                                "North Europe" = "#fac03d", 
#                                "Southern Europe" = "#db9b34", 
#                                "Western Europe" = "#c67915", 
#                                "South America" = "#ed6d3d",
#                                "North America" = "#e60012",
#                                "North Africa" = "#bba1cb", 
#                                "Western Africa" = "#ba79b1", 
#                                "Southern Africa" = "#8a1874", 
#                                "Eastern Africa" = "#422256"),
#                     na.value = "grey95") +
#   ggnewscale::new_scale_fill() + 
#   geom_fruit(
#     geom = geom_tile,
#     mapping = aes(y = tip, fill = plot_year),
#     width = 0.0025,
#     offset = 0.06
#   ) + 
#   scale_fill_manual(values = c("before 1990s" = "#eaf4fc", "1990s" = "#bccde9", "2000s" = "#8ea5d6", "2010s" = "#607dc3", "2020s" = "#0f2350"),
#                     na.value = "grey95")


lineage_nodes <- metadata_tree %>%
  filter(!is.na(group_2)) %>%
  group_by(group_2) %>%
  summarise(
    node = getMRCA(raxml_tree, tip),
    label1 = tip[1],  # 左端 tip
    label2 = tip[n()] # 右端 tip
  ) %>% 
  rename(lineage = group_2)

ggtree(ape::compute.brlen(raxml_tree, method = "Grafen"), size = 0.3, layout = "fan", open.angle = 30) %<+%
  (
    metadata_tree %>% 
      mutate(
        plot_year = case_when(
          collection_year < 1990 ~ 'before 1990s', 
          collection_year >= 1990 & collection_year < 2000 ~ '1990s', 
          collection_year >= 2000 & collection_year < 2010 ~ '2000s', 
          collection_year >= 2010 & collection_year < 2020 ~ '2010s', 
          collection_year >= 2020 ~ '2020s', 
          .default = NA
        ), 
        area = case_when(
          area == "Yunnan" ~ "East Asia", 
          .default = area
        )
      )
  ) + 
  geom_hilight(data = lineage_nodes, aes(node = node, fill = lineage), type = "roundrect") + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                                "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                                "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                     na.value = "transparent") + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = continent),
    width = 0.05,
    offset = 0.075
  ) + 
  scale_fill_manual(values = c("Asia" = "#8ECFC9", "Europe" = "#FFBE7A", "Africa" = "#bba1cb", "America" = "#FA7F6F"), na.value = "grey95") +
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = area),
    width = 0.05,
    offset = 0.1
  ) + 
  scale_fill_manual(values = c("East Asia" = "#4f6f46", 
                               "Western Asia" = "#99bcac", 
                               "Southeast Asia" = "#a9be7b", 
                               "South Asia" = "#c0d695", 
                               "North Asia" = "#b7d332", 
                               "Eastern Europe" = "#fce2c4", 
                               "Central Europe" = "#fedc5e", 
                               "North Europe" = "#fac03d", 
                               "Southern Europe" = "#db9b34", 
                               "Western Europe" = "#c67915", 
                               "South America" = "#ed6d3d",
                               "North America" = "#e60012",
                               "North Africa" = "#bba1cb", 
                               "Western Africa" = "#ba79b1", 
                               "Southern Africa" = "#8a1874", 
                               "Eastern Africa" = "#422256"),
                    na.value = "grey95") +
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = plot_year),
    width = 0.05,
    offset = 0.1
  ) + 
  scale_fill_manual(values = c("before 1990s" = "#eaf4fc", "1990s" = "#bccde9", "2000s" = "#8ea5d6", "2010s" = "#607dc3", "2020s" = "#0f2350"),
                    na.value = "grey95") + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = label_yunnan),
    width = 0.05,
    offset = 0.1
  ) + 
  scale_fill_manual(values = c("Yunnan" = "black"), na.value = "grey97") + 
  theme(
    legend.position = "none" # 去除图例
  )

# ggsave('./plot/tree/raxml_tree_col2_type_2.pdf', width = 8, height = 8)

