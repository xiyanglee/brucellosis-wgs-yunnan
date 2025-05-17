# BiocManager::install("treeio", lib = "/public/r_share_library/4.1")
library(tidyverse)
library(treeio)
library(ape)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')
load('./r_image/metadata_v3.RData')



# 读取共识树
beast_tree <- treeio::read.beast('./data/beast_gubbins_contig/run_1/beast_GTRIG_clean_core_filtered.tree')
# 读取采样树集合
beast_tree.sampled <- treeio::read.beast('./data/beast_gubbins_contig/run_1/beast_GTRIG_clean_core_subsampled2.trees')

# 保留指定 tip 子集
tips_to_keep <- metadata.phylo$assembly_accession[metadata.phylo$lineage_level_1 == '1']
beast_tree.l1 <- treeio::drop.tip(beast_tree, setdiff(beast_tree@phylo$tip.label, tips_to_keep))
beast_tree.sampled.l1 <- lapply(beast_tree.sampled, function(tree) {
  treeio::drop.tip(tree, setdiff(tree@phylo$tip.label, tips_to_keep))
})

# --------------------------------------------------
# 提取共识树节点 height 及分布
mrca.height <- beast_tree.l1[, c("node", "height", "height_0.95_HPD", "height_median", "height_range")]
mrca.height <- mrca.height %>% 
  mutate(
    height = as.numeric(height),
    # height_median = as.numeric(height_median),
    # height_HPD_lower = sapply(height_0.95_HPD, function(x) x[1]),
    # height_HPD_upper = sapply(height_0.95_HPD, function(x) x[2]),
    # height_range_lower = sapply(height_range, function(x) x[1]),
    # height_range_upper = sapply(height_range, function(x) x[2])
  ) %>%
  select(node, height) # , height_range_lower, height_HPD_lower, height_median, height_HPD_upper, height_range_upper)

# --------------------------------------------------
# 生成两两样本对 + 共识树 MRCA 编号
sample_ids <- beast_tree.l1@phylo$tip.label
sample_pairs <- combn(sample_ids, 2)

mrcas <- vapply(
  seq_len(ncol(sample_pairs)),
  function(i) getMRCA(beast_tree.l1@phylo, sample_pairs[, i]),
  FUN.VALUE = integer(1)
)

mrca.nodes <- data.frame(
  sample_1 = sample_pairs[1, ],
  sample_2 = sample_pairs[2, ],
  mrca_node = mrcas
)

rm(sample_pairs, mrcas)

# --------------------------------------------------
# 采样采 node.depth.edgelength 中 mrca.height$node 的 depth

# 定义函数：根据采样树，提取 mrca.height 中 node 对应的 MRCA 时间
MAX_DATE <- 2023.58333333333
get_mrca_heights_from_sampled_tree <- function(tree_phylo, node_vec) {
  all_depths <- node.depth.edgelength(tree_phylo)
  max_depth <- max(all_depths)
  # 提取目标节点
  selected_depths <- all_depths[node_vec]
  # max_depth <- max(selected_depths)
  # 转换成 "时间"（树顶到节点）
  height_times <- max_depth - selected_depths
  return(height_times)
}

# 批量处理所有采样树
sampled_heights_list <- lapply(beast_tree.sampled.l1, function(tree) {
  get_mrca_heights_from_sampled_tree(tree@phylo, mrca.height$node)
})

# 整合为矩阵
sample_matrix <- do.call(cbind, sampled_heights_list)
colnames(sample_matrix) <- paste0("height_", seq_len(ncol(sample_matrix)))

# 并到 mrca.height
mrca.height <- cbind(mrca.height, as.data.frame(sample_matrix))

rm(sampled_heights_list, sample_matrix)

# --------------------------------------------------
# 与 mrca.nodes 合并

mrca.nodes <- left_join(mrca.nodes, mrca.height, by = c('mrca_node' = 'node'))

# 合并地理/谱系信息
mrca.nodes <- mrca.nodes %>% 
  left_join(
    ., 
    metadata.phylo %>% 
      dplyr::select(assembly_accession, lineage_level_2, geographic_country, continent, geographic_province, collection_year) %>%
      rename_with(~paste0(., "_1"), .cols = c(lineage_level_2, geographic_country, continent, geographic_province, collection_year)), 
    by = c("sample_1" = "assembly_accession")
  ) %>% 
  left_join(
    ., 
    metadata.phylo %>% 
      dplyr::select(assembly_accession, lineage_level_2, geographic_country, continent, geographic_province, collection_year) %>%
      rename_with(~paste0(., "_2"), .cols = c(lineage_level_2, geographic_country, continent, geographic_province, collection_year)), 
    by = c("sample_2" = "assembly_accession")
  )

# save.image('./r_image/global_dynamic_trees_CI_base_info.RData')


# 设置 mrca_group、year_diff、geo_group
mrca.nodes.fit <- mrca.nodes %>%
  mutate(
    mrca_group = case_when(
      height >= 0  & height < 10   ~ "0-10 y",
      height >= 10  & height < 50  ~ "10-50 y",
      height >= 50  & height < 100 ~ "50-100 y",
      height >= 100 & height < 500 ~ "100-500 y",
      height >= 500               ~ ">500 y",
      TRUE ~ NA_character_
    ),
    year_diff = abs(collection_year_1 - collection_year_2),
    geo_group = case_when(
      continent_1 != continent_2 ~ "continental",
      geographic_country_1 != geographic_country_2 ~ "national",
      geographic_province_1 != geographic_province_2 ~ "regional",
      geographic_province_1 == geographic_province_2 ~ "local",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!(geo_group == "same_province" & geographic_province_1 == "" & geographic_province_2 == ""))

table(mrca.nodes.fit$mrca_group, mrca.nodes.fit$geo_group)

ref_group <- "regional"

rr_summary <- mrca.nodes.fit %>%
  group_by(mrca_group, geo_group) %>%
  summarise(pair_count = n(), .groups = "drop") %>%
  group_by(mrca_group) %>%
  mutate(
    total = sum(pair_count),
    prob = pair_count / total,
    ref_prob = prob[geo_group == ref_group],
    RR = prob / ref_prob
  )

ggplot(rr_summary, aes(x = RR, y = geo_group)) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(size = 3) +
  scale_x_log10() + 
  theme_classic(base_size = 14) + 
  coord_flip() + 
  facet_wrap(~ mrca_group, ncol = 1)


# 获取 height_1 到 height_100 的列名
height_cols <- paste0("height_", 1:100)

# 提取子数据并准备嵌套操作
mrca_long_list <- lapply(seq_along(height_cols), function(i) {
  col_i <- height_cols[i]
  
  mrca.nodes %>%
    # transmute( # same as mutate(.keep = "none")
    mutate(
      height = .data[[col_i]],
      geo_group = case_when(
        continent_1 != continent_2 ~ "continental",
        geographic_country_1  != geographic_country_2  ~ "national",
        geographic_province_1 != geographic_province_2 ~ "regional",
        geographic_province_1 == geographic_province_2 ~ "local",
        TRUE ~ NA_character_
      ),
      mrca_group = case_when(
        height >= 0  & height < 10   ~ "0-10 y",
        height >= 10  & height < 50  ~ "10-50 y",
        height >= 50  & height < 100 ~ "50-100 y",
        height >= 100 & height < 500 ~ "100-500 y", 
        height >= 500                ~ ">500 y",
        TRUE ~ NA_character_
      ), 
      geographic_province_1 = geographic_province_1, 
      geographic_province_2 = geographic_province_2, 
      .keep = "none" # 速度更快
    ) %>% 
    dplyr::filter(!(geo_group == "same_province" & geographic_province_1 == "" & geographic_province_2 == "")) %>% 
    dplyr::filter(!is.na(geo_group), !is.na(mrca_group)) %>% 
    group_by(mrca_group, geo_group) %>%
    summarise(pair_count = n(), .groups = "drop") %>%
    group_by(mrca_group) %>%
    mutate(
      total = sum(pair_count),
      prob = pair_count / total,
      ref_prob = prob[geo_group == ref_group],
      RR = prob / ref_prob
    ) %>%
    ungroup() %>%
    select(mrca_group, geo_group, RR)
})

# 合并 100 次模拟结果
rr_all <- bind_rows(mrca_long_list, .id = "replicate")  # replicate 1~100

# 计算每组的 Q05 和 Q95
rr_ci <- rr_all %>%
  group_by(mrca_group, geo_group) %>%
  summarise(
    RR_lower = quantile(RR, probs = 0.05, na.rm = TRUE),
    RR_upper = quantile(RR, probs = 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# 合并回 rr_summary
rr_summary <- rr_summary %>%
  left_join(rr_ci, by = c("mrca_group", "geo_group"))

rm(mrca_long_list, rr_all, rr_ci)

rr_summary <- rr_summary %>% 
  mutate(
    geo_group = factor(geo_group, levels = c("local", "regional", "national", "continental")), 
    mrca_group = factor(mrca_group, levels = c("0-10 y", "10-50 y", "50-100 y", "100-500 y", ">500 y"))
  )

ggplot(rr_summary, aes(x = RR, y = geo_group)) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +  # 参考线
  geom_errorbarh(aes(xmin = RR_lower, xmax = RR_upper), height = 0.1, color = "black", linewidth = 0.5) +  # 置信区间
  geom_point(size = 2.5, color = "black") +  # RR 点
  scale_x_log10() + 
  theme_classic(base_size = 14) + 
  coord_flip() + 
  facet_wrap(~ mrca_group, ncol = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave('./plot/geospatial_dynamic/lineage_1_RR_across_diff_spatial_scale_by_MRCA.pdf', width = 3, height = 9)

