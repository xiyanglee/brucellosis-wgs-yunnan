# BiocManager::install("treeio", lib = "/public/r_share_library/4.1")
library(tidyverse)
library(treeio)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



beast_tree <- treeio::read.beast('./data/beast_gubbins_contig/run_1/beast_GTRIG_clean_core_filtered.tree')
tips_to_keep  <- metadata.phylo$assembly_accession[metadata.phylo$lineage_level_1 == '1']
# beast_tree.l1 <- beast_tree
beast_tree.l1 <- treeio::drop.tip(beast_tree, setdiff(beast_tree@phylo$tip.label, tips_to_keep))

## 查看 beast 信息
get.fields(beast_tree.l1)
# "default.rate"          "default.rate_0.95_HPD" "default.rate_median"   "default.rate_range"    
# "height"                "height_0.95_HPD"       "height_median"         "height_range"          
# "length"                "length_0.95_HPD"       "length_median"         "length_range"         
# "posterior"             "state"                 "state.prob"            
# "state.rate"            "state.rate_0.95_HPD"   "state.rate_median"     "state.rate_range"      
# "state.set"             "state.set.prob"       


# ---------------------------
# 提取节点对应的 height 及分布
# 

## 测试：验证提取子树后 getMRCA 提取的 node 对应的 height 是否会变化
## 测试: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# getMRCA(beast_tree@phylo, c("202127", "202243"))
# # 1338
# data_mrca <- beast_tree[, c("node", "height")]
# data_mrca$height <- as.numeric(data_mrca$height)
# data_mrca$height[data_mrca$node == 1338]
# # [1] 48.97234
# 
# getMRCA(beast_tree.l1@phylo, c("202127", "202243"))
# # 1040
# data_mrca <- beast_tree.l1[, c("node", "height")]
# data_mrca$height <- as.numeric(data_mrca$height)
# data_mrca$height[data_mrca$node == 1040]
# # [1] 48.97234
## 测试: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

MAX_DATE <- 2023.58333333333
mrca.height <- beast_tree.l1[, c("node", "height", "height_0.95_HPD", "height_median", "height_range")]
mrca.height <- mrca.height %>% 
  mutate(
    height = as.numeric(height),
    height_median = as.numeric(height_median),
    height_HPD_lower = sapply(height_0.95_HPD, function(x) x[1]),
    height_HPD_upper = sapply(height_0.95_HPD, function(x) x[2]),
    height_range_lower = sapply(height_range, function(x) x[1]),
    height_range_upper = sapply(height_range, function(x) x[2])
  ) %>%
  select(node, height, height_range_lower, height_HPD_lower, height_median, height_HPD_upper, height_range_upper)


# ---------------------------
# 生成两两样本对，并获取样本对的 MRCA 编号
# 
sample_ids <- beast_tree.l1@phylo$tip.label
sample_pairs <- combn(sample_ids, 2)

# 获取 MRCA 编号
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


# ---------------------------
# 生成两两样本对，并获取样本对的 MRCA 编号
# 

## 实验: 给定分位数，返回符合分布的抽样结果
## 测试: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# val <- unname(unlist(data_mrca[1001, c(3:7)]))
# q <- c(0, 0.025, 0.5, 0.975, 1)
# # 拟合 spline CDF，并采样
# quantile_fun <- splinefun(q, val, method = "monoH.FC")
# samples <- quantile_fun(runif(100))
## 测试: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

## 注意，当 node 的 height 为 0 时
## 其 height_range_lower height_HPD_lower 等参数都为 NA
## 测试: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# p <- ggtree(beast_tree.l1, branch.length = "height") +
#   geom_tree()
# # 获取对应节点的坐标信息
# tree_data <- p$data
# node_target <- 395
# node_coords <- tree_data[tree_data$node == node_target, ]
# # 添加标记（用红色圆点 + 标签）
# p + 
#   geom_point(data = node_coords, aes(x = x, y = y), color = "red", size = 3) +
#   geom_text(data = node_coords, aes(x = x, y = y + 1, label = paste("Node", node)), 
#             color = "red", size = 3, fontface = "bold")
# rm(p, tree_data, node_target, node_coords)
## 测试: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 为 mrca.height 中每一行拟合并生成 100 个采样
generate_row_samples_from_specified_distribution <- function(row_vals, n = 100) {
  # 定义分位数位置
  q <- c(0, 0.025, 0.5, 0.975, 1)
  val <- unname(unlist(row_vals))
  
  # 如果节点的 height 为 0，各个参数都为 NA
  # 若 val 中含 NA 或长度不为 5，则返回 0 向量
  if (any(is.na(val)) || length(val) != 5) {
    return(rep(0, n))
  }
  
  quantile_fun <- splinefun(q, val, method = "monoH.FC")
  samples <- quantile_fun(runif(n))
  return(samples)
}

sample_matrix <- apply(
  mrca.height[, c("height_range_lower", "height_HPD_lower", "height_median", "height_HPD_upper", "height_range_upper")],
  1,
  generate_row_samples_from_specified_distribution,
  n = 100
)

# 转置为行对齐矩阵（行 = 节点，列 = height_1 到 height_100）
sample_matrix <- t(sample_matrix)
# 构造列名
colnames(sample_matrix) <- paste0("height_", seq_len(ncol(sample_matrix)))
# 合并到原始 mrca.height 中
mrca.height <- cbind(mrca.height, as.data.frame(sample_matrix))

rm(sample_matrix)


# save.image('./r_image/geospatial_dynamic_mrca_base_info.RData')

################################################################################

# load('./r_image/geospatial_dynamic_mrca_base_info.RData')

# ---------------------------
# 两两样本对合并地理、谱系信息
# 

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


mrca.nodes <- left_join(mrca.nodes, mrca.height, by = c('mrca_node' = 'node'))

## 将 MRCA 的最小值设为 0
# min(mrca.nodes$height); min(mrca.nodes$height_1); min(mrca.nodes$height_2)
# [1] 1.920581
# [1] 2.308308
# [1] 1.943951
# mrca.nodes <- mrca.nodes %>%
#   mutate(across(starts_with("height"), ~ .x - min(.x, na.rm = TRUE)))
# min(mrca.nodes$height); min(mrca.nodes$height_1); min(mrca.nodes$height_2)
# [1] 0
# [1] 0
# [1] 0

# # ---- Step B: 标记是否与谱系 1.3 相关 ----
# mrca.nodes <- mrca.nodes %>%
#   mutate(
#     involves_lineage = (lineage_level_2_1 == "1.3") | (lineage_level_2_2 == "1.3")
#   ) %>%
#   filter(involves_lineage)

mrca.nodes.fit <- mrca.nodes %>%
  mutate(
    mrca_group = case_when(
      height >= 0  & height < 10   ~ "0-10 y",
      height >= 10  & height < 50  ~ "10-50 y",
      height >= 50  & height < 100 ~ "50-100 y",
      height >= 100 & height < 500 ~ "100-500 y",
      height >= 100                ~ ">500 y",
      TRUE                        ~ NA_character_
    ),
    year_diff = abs(collection_year_1 - collection_year_2),
    geo_group = case_when(
      continent_1 != continent_2 ~ "continental",
      geographic_country_1  != geographic_country_2  ~ "national",
      geographic_province_1 != geographic_province_2 ~ "regional",
      geographic_province_1 == geographic_province_2 ~ "local",
      TRUE ~ NA_character_
    )
  ) %>%
  # filter(!is.na(geo_group) & year_diff <= 5)
  dplyr::filter(!(geo_group == "same_province" & geographic_province_1 == "" & geographic_province_2 == ""))


# set.seed(12345)  # 确保结果可复现（随机分配）
# mrca.nodes.fit <- mrca.nodes %>%
#   mutate(
#     mrca_group = case_when(
#       height >= 0  & height < 10   ~ "0-10 y",
#       height >= 10  & height < 50  ~ "10-50 y",
#       height >= 50  & height < 100 ~ "50-100 y",
#       height >= 100 & height < 500 ~ "100-500 y", 
#       height >= 500               ~ ">500 y",
#       TRUE                        ~ NA_character_
#     ), 
#     year_diff = abs(collection_year_1 - collection_year_2),
#     
#     # 初始化 geo_group
#     geo_group = case_when(
#       continent_1 != continent_2 ~ "inter_continent",
#       geographic_country_1  != geographic_country_2  ~ "diff_country",
#       geographic_province_1 != geographic_province_2 ~ "intra_country_diff_province",
#       geographic_province_1 == geographic_province_2 ~ "same_province",
#       TRUE ~ NA_character_
#     ),
#     
#     # 补充特殊情况：两个省份均为空字符串
#     geo_group = if_else(
#       geographic_province_1 == "" & geographic_province_2 == "",
#       if_else(runif(n()) < 0.5, "same_province", "intra_country_diff_province"),
#       geo_group
#     )
#   )

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

