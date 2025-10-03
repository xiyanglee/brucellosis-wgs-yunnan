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

raxml_tree <- ape::read.tree('./data/raxml_gubbins_contig/RAxML_bestTree.raxml_tree')
# 重新设定外类群
raxml_tree <- root(raxml_tree, outgroup = 'GCA_000369945.1', resolve.root = TRUE)

plot(raxml_tree, show.node.label = T)
ggtree(raxml_tree)

fasta_sparse_data <- import_fasta_sparse_nt('./data/raxml_gubbins_contig/filtered_branch.core.aln')
keep_tip <- !(colnames(fasta_sparse_data$snp.matrix) %in% drop_tip) & (colnames(fasta_sparse_data$snp.matrix) %in% raxml_tree$tip.label)
fasta_sparse_data$snp.matrix <- fasta_sparse_data$snp.matrix[, keep_tip]

fasta_sparse_data_1 <- optimise_prior(fasta_sparse_data, type = "symmetric", hc.method = "ward", grid.interval = c(5e-9, 1))
raxml_hclust <- phytools::midpoint.root(raxml_tree)
best_partition_1 <- best_baps_partition(fasta_sparse_data_1, raxml_hclust)

group_data_1 <- as.data.frame(best_partition_1) %>% 
  rownames_to_column()
colnames(group_data_1) <- c('tip', 'group_1')
group_data_1$group_1 <- as.factor(group_data_1$group_1)

ggtree(raxml_tree) %<+% group_data_1 + geom_tippoint(aes(color = group_1), size = 2) + 
  ggsci::scale_color_bmj()

group_data_1 <- group_data_1 %>% 
  mutate(group_1 = case_when(
    group_1 == 2 ~ '1', 
    group_1 == 3 ~ '2', 
    group_1 == 1 ~ '3'
  ))

fasta_sparse_data_2 <- optimise_prior(fasta_sparse_data, type = "baps", hc.method = "ward", grid.interval = c(5e-9, 1), n.cores = 8)
# [1] "Optimised hyperparameter: 0.039"
best_partition_2 <- best_baps_partition(fasta_sparse_data_2, raxml_hclust)

group_data_2 <- as.data.frame(best_partition_2) %>% 
  rownames_to_column()
colnames(group_data_2) <- c('tip', 'group_4')
group_data_2$group_4 <- as.factor(group_data_2$group_4)

ggtree(raxml_tree) %<+% 
  group_data_2 + 
  geom_tippoint(aes(color = group_4), size = 2) + 
  ggsci::scale_color_npg()


boot.result <- boot_fast_baps(fasta_sparse_data_2, n.cores = 8, quiet = FALSE)
dendro <- as.dendrogram(fast_baps(fasta_sparse_data_2))
gplots::heatmap.2(boot.result, dendro, dendro, tracecol = NA)



# 转换 bootstrap 矩阵为长表格
boot_long <- as.data.frame(as.table(as.matrix(boot.result)))
colnames(boot_long) <- c("Sample1", "Sample2", "BootstrapValue")

# 添加分组信息
boot_long <- boot_long %>%
  left_join(group_data_2, by = c("Sample1" = "tip")) %>%
  rename(Group1 = group_4) %>%
  left_join(group_data_2, by = c("Sample2" = "tip")) %>%
  rename(Group2 = group_4)

# 计算每组对的平均值
group_summary <- boot_long %>%
  # filter(Sample1 != Sample2) %>%  # 去掉自比较
  group_by(Group1, Group2) %>%
  summarise(AverageBootstrap = mean(BootstrapValue, na.rm = TRUE), .groups = "drop")

# 转换为矩阵形式
group_matrix <- reshape2::dcast(group_summary, Group1 ~ Group2, value.var = "AverageBootstrap")
rownames(group_matrix) <- group_matrix$Group1
group_matrix <- group_matrix[, -1]

ggplot(group_summary, aes(x = Group1, y = Group2, fill = AverageBootstrap)) +
  geom_tile(color = "white") +  # 每个格子用白色边框分隔
  viridis::scale_fill_viridis(option = "magma", name = "Avg Bootstrap") + 
  labs(title = "Group-wise Average Bootstrap Heatmap",
       x = "Group 1",
       y = "Group 2",
       fill = "Avg Bootstrap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

group_summary %>%
  filter(Group1 != Group2 & AverageBootstrap > 50)
#   Group1 Group2 AverageBootstrap
#   <fct>  <fct>             <dbl>
# 1 2      3                  52  
# 2 3      2                  52  
# 3 4      7                 100  
# 4 5      10                 84.0
# 5 7      4                 100  
# 6 10     5                  84.0





ggtree(raxml_tree) %<+% 
  (
    group_data_2 %>% 
      mutate(mark = case_when(
        group_4 == 16 ~ 'mark', 
        .default = NA
      ))
  ) + 
  geom_tippoint(aes(color = mark), size = 2) + 
  scale_color_manual(values = c("mark" = "gold"), na.value = "transparent")

group_data_2 <- group_data_2 %>% 
  mutate(group_2 = case_when(
    group_4 == 5 ~ '1.1', 
    group_4 == 6 ~ '1.2', 
    group_4 %in% c(4, 7, 10) ~ '1.3', 
    group_4 == 8 ~ '2.1', 
    group_4 == 1 ~ '2.2', 
    group_4 == 9 ~ '2.3', 
    group_4 == 3 ~ '3.1', 
    group_4 == 2 ~ '3.2', 
    .default = NA
  ), 
  group_1  = case_when(
    group_2 %in% c('1.1', '1.2', '1.3') ~ '1', 
    group_2 %in% c('2.1', '2.2', '2.3') ~ '2', 
    group_2 %in% c('3.1', '3.2') ~ '3'
  ))



# 转换 bootstrap 矩阵为长表格
boot_long_2 <- as.data.frame(as.table(as.matrix(boot.result)))
colnames(boot_long_2) <- c("Sample1", "Sample2", "BootstrapValue")

# 添加分组信息
boot_long_2 <- boot_long_2 %>%
  left_join(dplyr::select(group_data_2, tip, group_2), by = c("Sample1" = "tip")) %>%
  rename(Group1 = group_2) %>%
  left_join(dplyr::select(group_data_2, tip, group_2), by = c("Sample2" = "tip")) %>%
  rename(Group2 = group_2)

# 计算每组对的平均值
group_summary_2 <- boot_long_2 %>%
  # filter(Sample1 != Sample2) %>%  # 去掉自比较
  group_by(Group1, Group2) %>%
  summarise(AverageBootstrap = mean(BootstrapValue, na.rm = TRUE), .groups = "drop")

ggplot(group_summary_2, aes(x = Group1, y = Group2, fill = AverageBootstrap)) +
  geom_tile(color = "white") +  # 每个格子用白色边框分隔
  geom_text(aes(label = ifelse(AverageBootstrap > 50, round(AverageBootstrap, 1), "")),
            color = "white", size = 3) +  # 显示大于 50 的值
  # viridis::scale_fill_viridis(option = "rocket", name = "Avg Bootstrap", direction = -1) + 
  scale_fill_gradientn(colors = c("grey95", "#a6bddb", "#034e7b"), name = "Avg. bootstrap") +  # 自定义配色
  labs(title = NULL, # "Lineage-wise Average Bootstrap Heatmap"
       x = "Lineage 1",
       y = "Lineage 2",
       fill = "Avg. bootstrap") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  coord_fixed()

# ggsave('./plot/tree/baps_bootstrap_heatmap.pdf', width = 5, height = 4)



ggtree(raxml_tree) %<+% group_data_2 + geom_tippoint(aes(color = group_2), size = 2) + 
  scale_color_manual(values = c("1.1" = "#FBCBD2", "1.2" = "#ee7959", "1.3" = "#ab1d22", 
                                "2.1" = "#aed0ee", "2.2" = "#2ca9e1", "2.3" = "#2e59a7", 
                                "3.1" = "#c0d695", "3.2" = "#4f6f46"), 
                     na.value = "transparent")

# save.image('./r_image/cluster_tree_full.RData')

rm(list = setdiff(ls(), c("group_data_2", "raxml_tree")))

# save.image('./r_image/cluster_tree.RData')


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

ggtree(raxml_tree, linewidth = 0.25) %<+%
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
  geom_tiplab(aes(color = label_yunnan), as_ylab = TRUE, align = T, linesize = 1, linetype = 1, offset = 0.001) + 
  scale_color_manual(values = c("Yunnan" = "grey92"), na.value = "transparent") + 
  theme(
    axis.text.y = element_blank()  # 移除 Y 轴刻度文本
  ) +
  ggnewscale::new_scale_color() + 
  geom_tippoint(aes(color = group_2), size = 0.5) + 
  scale_color_manual(values = c("1.1" = "#FBCBD2", "1.2" = "#ee827c", "1.3" = "#c9171e", 
                                "2.1" = "#aed0ee", "2.2" = "#2ca9e1", "2.3" = "#2e59a7", 
                                "3.1" = "#c0d695", "3.2" = "#4f6f46"), 
                     na.value = "transparent") + 
  ggnewscale::new_scale_color() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = continent),
    width = 0.0025,
    offset = 0.05
  ) + 
  scale_fill_manual(values = c("Asia" = "#8ECFC9", "Europe" = "#FFBE7A", "Africa" = "#bba1cb", "America" = "#FA7F6F"),
                    na.value = "grey95") +
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = area),
    width = 0.0025,
    offset = 0.06
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
    width = 0.0025,
    offset = 0.06
  ) + 
  scale_fill_manual(values = c("before 1990s" = "#eaf4fc", "1990s" = "#bccde9", "2000s" = "#8ea5d6", "2010s" = "#607dc3", "2020s" = "#0f2350"),
                    na.value = "grey95") +
  NULL

# ggsave('./plot/tree/raxml_tree.pdf', width = 8, height = 8)



metadata_tree %>%
  group_by(group_1, continent) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group_1) %>%
  mutate(proportion = count / sum(count)) %>% 
  ggplot(aes(x = group_1, y = proportion, fill = continent)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Lineage",
    y = "Proportion",
    fill = "Continent",
    title = NULL
  ) +
  theme_classic() + 
  scale_fill_manual(values = c("Asia" = "#8ECFC9", "Europe" = "#FFBE7A", "Africa" = "#bba1cb", "America" = "#FA7F6F"),
                    na.value = "grey95")

metadata_tree %>% 
  dplyr::select(group_2, continent) %>% 
  na.omit() %>% 
  group_by(group_2, continent) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group_2) %>%
  mutate(proportion = count / sum(count)) %>% 
  ggplot(aes(x = group_2, y = proportion, fill = continent)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    x = "Lineage",
    y = "Proportion",
    fill = "Continent",
    title = NULL
  ) +
  theme_classic() + 
  scale_fill_manual(values = c("Asia" = "#8ECFC9", "Europe" = "#FFBE7A", "Africa" = "#bba1cb", "America" = "#FA7F6F"),
                    na.value = "grey95")

# ggsave('./plot/tree/continent_proportion_in_lineage.pdf', width = 4, height = 3.5)

metadata_tree %>% 
  mutate(
    area = case_when(
      area == "Yunnan" ~ "Yunnan (This study)", # East Asia
      .default = area
    )
  ) %>%
  dplyr::select(group_2, area) %>% 
  na.omit() %>% 
  group_by(group_2, area) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group_2) %>% 
  mutate(proportion = count / sum(count)) %>% 
  ggplot(aes(x = group_2, y = proportion, fill = area)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    x = "Lineage",
    y = "Area",
    fill = "Continent",
    title = NULL
  ) +
  theme_classic() + 
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
                               "Eastern Africa" = "#422256", 
                               "Yunnan (This study)" = "darkgreen"),
                    na.value = "grey95")

# ggsave('./plot/tree/area_proportion_in_lineage.pdf', width = 4, height = 3.5)
