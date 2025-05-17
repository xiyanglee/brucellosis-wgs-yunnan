# BiocManager::install("treeio", lib = "/public/r_share_library/4.1")
library(tidyverse)

setwd('/Users/xiyangli/Lab/Project/Brucella_WGS_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_data/BactDating.RData')



time_tree <- dated_tree$tree
# 获取要保留的 tip
sample_ids  <- metadata.phylo$assembly_accession[metadata.phylo$lineage_level_1 == '1' & !is.na(metadata.phylo$continent)]
sample_ids  <- intersect(time_tree$tip.label, sample_ids)


# ---------------------------
# 生成两两样本对，并获取样本对的 MRCA 编号
# 
sample_pairs <- combn(sample_ids, 2)

# 获取 MRCA 编号
mrcas <- vapply(
  seq_len(ncol(sample_pairs)),
  function(i) as.integer(getMRCA(time_tree, sample_pairs[, i])),
  FUN.VALUE = integer(1)
)

mrca.nodes <- data.frame(
  sample_1 = sample_pairs[1, ],
  sample_2 = sample_pairs[2, ],
  mrca_node = mrcas
)

rm(sample_pairs, mrcas)

# ---------------------------
# 获取样本对的 MRCA 编号
# 
mrca.nodes$date <- allDates(time_tree)[mrca.nodes$mrca_node]

time_tree.samples <- extractSample(dated_tree, 100)
for (i in seq_along(time_tree.samples)) {
  mrca.nodes[[paste0("date_", i)]] <- allDates(time_tree.samples[[i]])[mrca.nodes$mrca_node]
}; rm(i)

# plot(mrca.nodes$height, mrca.nodes$height_2)

mrca.nodes %>%
  select(starts_with("date_")) %>%
  unlist(use.names = FALSE) %>%
  max(na.rm = TRUE)
# [1] 2023.694

max(metadata.phylo$collection_date_float, na.rm = TRUE)
# 2023.583

max(mrca.nodes$date)
# [1] 2015.761

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

mrca.nodes.fit <- mrca.nodes %>%
  mutate(
    mrca_group = case_when(
      date >= 2005 ~ "0-10 y",
      # date >= 1965 ~ "10-50 y",
      date >= 1915 ~ "10-100 y",
      date >= 1515 ~ "100-500 y",
      date < 1515  ~ ">500 y",
      TRUE         ~ NA_character_
    ),
    geo_group = case_when(
      continent_1 != continent_2 ~ "continental",
      geographic_country_1  != geographic_country_2  ~ "national",
      geographic_province_1 != geographic_province_2 ~ "regional",
      geographic_province_1 == geographic_province_2 ~ "local",
      TRUE ~ NA_character_
    ), 
    year_diff = abs(collection_year_1 - collection_year_2),
  ) %>%
  filter(!is.na(geo_group) & year_diff <= 2) %>%
  dplyr::filter(!(geo_group == "local" & geographic_province_1 == "" & geographic_province_2 == "")) %>% 
  dplyr::filter(!(geo_group == "regional" & geographic_country_1 == "" & geographic_country_2 == "")) %>% 
  dplyr::filter(!(geo_group == "national" & continent_1 == "" & continent_2 == ""))

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
date_cols <- paste0("date_", 1:100)

# 提取子数据并准备嵌套操作
mrca_long_list <- lapply(seq_along(date_cols), function(i) {
  col_i <- date_cols[i]
  
  mrca.nodes %>%
    # transmute( # same as mutate(.keep = "none")
    mutate(
      date = .data[[col_i]],
      geo_group = case_when(
        continent_1 != continent_2 ~ "continental",
        geographic_country_1  != geographic_country_2  ~ "national",
        geographic_province_1 != geographic_province_2 ~ "regional",
        geographic_province_1 == geographic_province_2 ~ "local",
        TRUE ~ NA_character_
      ),
      mrca_group = case_when(
        date >= 2005 ~ "0-10 y",
        # date >= 1965 ~ "10-50 y",
        date >= 1915 ~ "10-100 y",
        date >= 1515 ~ "100-500 y",
        date < 1515  ~ ">500 y",
        TRUE         ~ NA_character_
      ), 
      year_diff = abs(collection_year_1 - collection_year_2),
      geographic_province_1 = geographic_province_1, 
      geographic_province_2 = geographic_province_2, 
      geographic_country_1 = geographic_country_1, 
      geographic_country_2 = geographic_country_2, 
      continent_1 = continent_1, 
      continent_2 = continent_2, 
      .keep = "none" # 速度更快
    ) %>% 
    filter(!is.na(geo_group) & year_diff <= 2) %>% 
    dplyr::filter(!(geo_group == "local" & geographic_province_1 == "" & geographic_province_2 == "")) %>% 
    dplyr::filter(!(geo_group == "regional" & geographic_country_1 == "" & geographic_country_2 == "")) %>% 
    dplyr::filter(!(geo_group == "national" & continent_1 == "" & continent_2 == "")) %>% 
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
    mrca_group = factor(mrca_group, levels = c("0-10 y", "10-100 y", "100-500 y", ">500 y"))
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

