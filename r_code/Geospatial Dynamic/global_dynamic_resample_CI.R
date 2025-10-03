library(tidyverse)
library(treeio)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

# 加载数据
load('./r_image/metadata_v3.RData')

beast_tree <- treeio::read.beast('./data/beast_gubbins_contig/run_1/beast_GTRIG_clean_core_filtered.tree')
tips_to_keep  <- metadata.phylo$assembly_accession[metadata.phylo$lineage_level_1 == '1']
beast_tree.l1 <- treeio::drop.tip(beast_tree, setdiff(beast_tree@phylo$tip.label, tips_to_keep))

# 取出MRCA节点height
MAX_DATE <- 2023.58333333333
mrca.height <- beast_tree.l1[, c("node", "height")]
mrca.height <- mrca.height %>% mutate(height = as.numeric(height))

# 两两组合
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

# 合并地理、年代信息
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
  ) %>%
  left_join(mrca.height, by = c('mrca_node' = 'node'))

# 加上分组
mrca.nodes.fit <- mrca.nodes %>%
  mutate(
    mrca_group = case_when(
      height >= 0  & height < 10   ~ "0-10 y",
      height >= 10  & height < 50  ~ "10-50 y",
      height >= 50  & height < 100 ~ "50-100 y",
      # height >= 100 & height < 500 ~ "100-500 y",
      height >= 100               ~ ">100 y",
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
  dplyr::filter(!(geo_group == "same_province" & geographic_province_1 == "" & geographic_province_2 == "")) %>% 
  dplyr::filter(!is.na(geo_group) & !is.na(mrca_group) & year_diff <= 5)


table(mrca.nodes.fit$mrca_group, mrca.nodes.fit$geo_group)
#          continental local national regional
# >100 y          9613 10091    53276    14076
# 0-10 y             0   169        0        6
# 10-50 y            4  3158       76      925
# 50-100 y         208   283      253      550

# ---------------------------
# 计算原始 RR 点估计
# ---------------------------
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
  ) %>%
  ungroup()

# ---------------------------
# Bootstrap 置信区间计算
# ---------------------------

set.seed(1234)
n_bootstrap <- 1000

# 所有样本
all_samples <- unique(c(mrca.nodes$sample_1, mrca.nodes$sample_2))

# 保存每次bootstrap的RR
bootstrap_results <- vector("list", length = n_bootstrap)

# 设定需要保证组合完整的水平
all_mrca_groups <- unique(mrca.nodes.fit$mrca_group)
all_geo_groups  <- unique(mrca.nodes.fit$geo_group)

for (b in seq_len(n_bootstrap)) {
  # 有放回bootstrap采样样本ID
  resampled_samples <- sample(all_samples, size = length(all_samples), replace = TRUE)
  resampled_samples <- unique(resampled_samples)
  
  # 筛选pair中两个样本都在重采样集合内的
  mrca.boot <- mrca.nodes.fit %>%
    filter(sample_1 %in% resampled_samples & sample_2 %in% resampled_samples)
  
  if (nrow(mrca.boot) == 0) next
  
  boot_rr <- mrca.boot %>% 
    ## 平衡每个 mrca_group 的数量
    ## method 1: 对所有超过数量的组进行抽样
    # group_by(mrca_group) %>%
    # group_modify(~ {
    #   n_group <- nrow(.x)
    #   n_mean <- mean(table(mrca.boot$mrca_group))  # 所有mrca_group的样本数均值
    #   n_mean <- 500 # 手动设置抽样阈值
    #   if (n_group > n_mean) {
    #     .x <- slice_sample(.x, n = round(n_mean))  # 超过均值则随机抽样
    #   }
    #   .x
    # }) %>%
    # ungroup() %>%
    ## method 2: 单独对指定组进行抽样
    group_by(mrca_group) %>%
    group_modify(~ {
      if (.y$mrca_group == ">100 y") {
        # 只对 ">500 y"组进行抽样
        slice_sample(.x, n = 100)
      } else if (.y$mrca_group == "50-100 y") {
        slice_sample(.x, n = 100)
      } else if (.y$mrca_group == "0-10 y") {
        slice_sample(.x, n = 100)
      } else {
        # 其他组不变
        .x
      }
    }) %>%
    ungroup() %>%
    ## ----------------------------------
    group_by(mrca_group, geo_group) %>%
    summarise(pair_count = n(), .groups = "drop") %>% 
    # 补全所有组合，缺失的设为0
    complete(mrca_group = all_mrca_groups, geo_group = all_geo_groups, fill = list(pair_count = 0)) %>% 
    # 把 pair_count = 0 的地方设为 1
    mutate(pair_count = ifelse(pair_count == 0, 1, pair_count)) %>% 
    # ----------------------------------
    group_by(mrca_group) %>%
    mutate(
      total = sum(pair_count),
      prob = pair_count / total,
      ref_prob = prob[geo_group == ref_group],
      RR = prob / ref_prob
    ) %>%
    ungroup() %>%
    select(mrca_group, geo_group, RR)
  
  bootstrap_results[[b]] <- boot_rr
}

# 合并所有bootstrap replicate
rr_boot_all <- bind_rows(bootstrap_results, .id = "replicate")

# 计算2.5%和97.5%的置信区间
rr_ci <- rr_boot_all %>%
  group_by(mrca_group, geo_group) %>%
  summarise(
    RR_lower = quantile(RR, probs = 0.025, na.rm = TRUE),
    RR_upper = quantile(RR, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# 合并回原rr_summary
rr_summary <- rr_summary %>%
  left_join(rr_ci, by = c("mrca_group", "geo_group"))

# ---------------------------
# 绘制最终结果
# ---------------------------

rr_summary <- rr_summary %>% 
  mutate(
    geo_group = factor(geo_group, levels = c("local", "regional", "national", "continental")),
    mrca_group = factor(mrca_group, levels = c("0-10 y", "10-50 y", "50-100 y", ">100 y"))
  )

# save.image('./r_image/global_dynamic_resample_CI.RData')

# load('./r_image/global_dynamic_resample_CI.RData')

ggplot(rr_summary, aes(x = RR, y = geo_group)) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = RR_lower, xmax = RR_upper), height = 0.1, color = "black", linewidth = 1) +
  geom_point(size = 2.8, color = "black") +
  scale_x_log10() + 
  # annotation_logticks(sides = "l") + 
  theme_classic(base_size = 14) +
  coord_flip() +
  facet_wrap(~ mrca_group, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 基础数据处理
rr_plot <- rr_summary %>%
  mutate(
    geo_group = factor(geo_group, levels = c("local", "regional", "national", "continental")),
    mrca_group = factor(mrca_group, levels = c("0-10 y", "10-50 y", "50-100 y", ">100 y")),
    bg_color = case_when(
      geo_group == "regional" ~ "#e6e4e3",
      geo_group == "local" ~ "#edf5e6",
      geo_group %in% c("national", "continental") ~ "#d2d1eb",
      TRUE ~ NA_character_
    ),
    point_shape = case_when(
      geo_group == "regional" ~ 0, # 0  空心正方形 
      geo_group == "local" ~ 17,   # 15 实心正方形 
      .default = 17                # 17 实心三角形
    )
  )

# 补充完整组合
full_combinations <- expand.grid(
  mrca_group = levels(rr_plot$mrca_group),
  geo_group = levels(rr_plot$geo_group)
)

# 合并，打上是否为新增缺失补齐行的标签
rr_plot_full <- full_combinations %>%
  left_join(rr_plot, by = c("mrca_group", "geo_group")) %>%
  mutate(
    is_missing = if_else(is.na(RR), TRUE, FALSE),
    bg_color = case_when(
      geo_group == "regional" ~ "#e6e4e3",
      geo_group == "local" ~ "#edf5e6",
      geo_group %in% c("national", "continental") ~ "#d2d1eb",
      TRUE ~ NA_character_
    ),
    point_shape = case_when(
      is_missing & geo_group %in% c("national", "continental") ~ 2,  # 空心三角
      is_missing & geo_group == "regional" ~ 0,
      is_missing & geo_group == "local" ~ 15,
      TRUE ~ point_shape
    ),
    RR = if_else(is_missing, 0.003, RR),
    # 保持已有RR_lower/RR_upper，新增补充行才设为NA
    RR_lower = if_else(is_missing, NA_real_, RR_lower),
    RR_upper = if_else(is_missing, NA_real_, RR_upper)
  )

# 正式绘图
ggplot(rr_plot_full, aes(x = RR, y = geo_group)) +
  
  # 背景色
  geom_tile(
    aes(x = 1, y = geo_group, fill = bg_color),
    width = Inf, height = 1,
    inherit.aes = FALSE, alpha = 0.4
  ) +
  
  # 参考线
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  
  # 置信区间，只画 is_missing==FALSE 的
  geom_errorbarh(
    data = rr_plot_full %>% filter(!is_missing, geo_group != "regional"),
    aes(xmin = RR_lower, xmax = RR_upper),
    height = 0.15, color = "black", linewidth = 0.5
  ) +
  
  # 点，按shape控制
  geom_point(aes(shape = factor(point_shape)), size = 2.5, color = "black", fill = "white", stroke = 0.75) +
  scale_shape_manual(values = c(`0` = 0, `15` = 15, `17` = 17, `2` = 2), guide = "none") +
  
  coord_flip() + 
  
  # x轴log10比例 + 正常标签
  scale_x_log10(
    limits = c(0.002, 500),
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100"),
    expand = c(0, 0)
  ) + 
  scale_y_discrete(expand = c(0, 0)) + 
  scale_fill_identity() +
  
  facet_wrap(~ mrca_group, ncol = 1) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
    # strip.text = element_text(face = "bold")
  ) + 
  labs(
    x = "Relative Risk",
    y = NULL
  )

# ggsave('./plot/geospatial_dynamic/lineage_1_global_RR_bootstrap_CI_plot.pdf', width = 3, height = 9)

# write.csv(rr_summary, "./result/geospatial_dynamic/lineage_1_global_RR_bootstrap_CI_summzry.csv", row.names = FALSE)

write.csv(rr_plot_full, "./r_plot_data/figure_5a_rr_plot_full.csv", row.names = FALSE)
