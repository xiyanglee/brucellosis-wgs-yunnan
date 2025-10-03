# BiocManager::install("treeio", lib = "/public/r_share_library/4.1")
# 用于 st_sample() 在计算随机点时，坐标范围沿着大圆计算（即没有考虑球面上的经纬度）
# install.packages("lwgeom", lib = "/public/r_share_library/4.1")
library(tidyverse)
library(treeio)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



geo_sf.province <- read_sf('../../metadata/中国_省.geojson') %>% st_make_valid() %>% dplyr::filter(name != "境界线")
geo_sf.city     <- read_sf('../../metadata/中国_市.geojson') %>% st_make_valid() %>% dplyr::filter(name != "境界线")
geo_sf.county   <- read_sf('../../metadata/中国_县.geojson') %>% st_make_valid() %>% dplyr::filter(name != "境界线")

geo_sf <- bind_rows(geo_sf.province, geo_sf.city, geo_sf.county) %>% distinct()
rm(geo_sf.province, geo_sf.city, geo_sf.county)

head(sort(table(geo_sf$name), decreasing = TRUE), 32)
# 境界线 市中区 鼓楼区   城区 新华区   郊区 铁西区   东区 南山区 南沙区 向阳区 
#     15      4      4      3      3      3      3      2      2      2      2 
# 和平区 城中区 城关区 宝山区 新城区 普陀区 朝阳区 桥西区 永定区 江北区 河东区 
#      2      2      2      2      2      2      2      2      2      2      2 
# 海州区 白云区 西安区 西湖区 通州区 铁东区 长安区 青山区 龙华区 丁青县 
#      2      2      2      2      2      2      2      2      2      1 

geo_sf <- geo_sf %>%
  mutate(
    centroid = st_centroid(geometry), # 计算质心
    lng = st_coordinates(centroid)[, 1], # 提取经度
    lat = st_coordinates(centroid)[, 2]  # 提取纬度
  ) %>%
  select(-c(centroid, gb)) %>%
  st_drop_geometry()

metadata.china <- metadata.raw.genebank %>% 
  dplyr::filter(geographic_country == "China")

sum(is.na(metadata.china$geographic_province))
# [1] 8

metadata.china <- metadata.china[!is.na(metadata.china$geographic_province), ]

metadata.china <- metadata.china %>% 
  mutate(geographic_province_cn = case_match(
    geographic_province, 
    "Kezuo Central Banner of Tongliao City, Tongliao, Inner Mongolia Autonomous region" ~ "科尔沁左翼中旗", 
    "Xinjiang"                    ~ "新疆维吾尔自治区", 
    "Hebei"                       ~ "河北省", 
    "Hohhot"                      ~ "呼和浩特市", 
    "Hulunbuir"                   ~ "呼伦贝尔市", 
    "Beijing"                     ~ "北京市", 
    "Qinghai"                     ~ "青海省", 
    "Tibetan Autonomous County of Tianzhu"                           ~ "天祝藏族自治县", 
    "Baotou City, Inner Mongolia Autonomous Region"                  ~ "包头市",                  
    "Ulanqab"                     ~ "乌兰察布市", 
    "Jarud Banner of Tongliao City, Inner Mongola Autonomous region" ~ "扎鲁特旗", 
    "Kailu Count of Tongliao City, Inner Mongolia Autonomous Region" ~ "开鲁县", 
    "hohhot"                      ~ "呼和浩特市", 
    "Inner Mongolia"              ~ "内蒙古自治区", 
    "Harbin, Northeastern China"  ~ "哈尔滨市", 
    "Hainan"                      ~ "海南省", 
    "Qinghai Province"            ~ "青海省", 
    "Inner Mongoalia Autonoumous Region"                             ~ "内蒙古自治区", 
    "Yangzhou"                    ~ "扬州市", 
    "Henan"                       ~ "河南省", 
    "Zibo"                        ~ "淄博市", 
    "Heilongjiang"                ~ "黑龙江省", 
    "Jilin"                       ~ "吉林省", 
    "Inner Mongolia, Ulanqab"     ~ "乌兰察布市", 
    "HaiNan"                      ~ "海南省",                                                      
    "LiaoNing"                    ~ "辽宁省",                                                      
    "JiaoZhou"                    ~ "胶州市", 
    "NeiMonggal"                  ~ "内蒙古自治区", 
    "Shandong"                    ~ "山东省",                                                      
    "XinJiang"                    ~ "新疆维吾尔自治区",                                                      
    "Jilin Province"              ~ "吉林省",                                                      
    "Shihezi Xinjiang"            ~ "石河子市",                                                      
    "Gansu"                       ~ "甘肃省", 
    .default = NA
  ))

metadata.china <- left_join(metadata.china, geo_sf, by = c("geographic_province_cn" = "name"))

# 取两个表格列的交集并合并
common_cols    <- intersect(colnames(metadata.china), colnames(metadata.sample))
metadata.china <- bind_rows(
  metadata.china[, common_cols],
  metadata.sample[, common_cols]
)
rm(common_cols)

# ------------------------------
# 计算测序样本两两之间的各变量差值
# ------------------------------

# 生成所有不重复的 i, j 组合
# combn(seq_len(n), 2) 返回一个 2 x (n*(n-1)/2) 的矩阵
# 每一列对应一对样本索引 (i, j)，确保了不重复且 i < j
# t() 将其转置后得到 n/2 行、2 列，更直观
n         <- nrow(metadata.china)
pairs_mat <- t(combn(seq_len(n), 2))  # 变为 nC2 行、2 列

# 从 pairs_mat 中分离 i, j
i_idx <- pairs_mat[, 1]
j_idx <- pairs_mat[, 2]

## 逐列向量化计算距离
## distHaversine() 支持将坐标输入为两列矩阵（每行对应一组坐标），会返回一个向量，向量长度与行数（即对数）匹配
dist_m <- distHaversine(
  cbind(metadata.china$lng[i_idx], metadata.china$lat[i_idx]),
  cbind(metadata.china$lng[j_idx], metadata.china$lat[j_idx])
)
dist_km <- dist_m / 1000

metadata.pairwise <- data.frame(
  tip_pair  = paste0(metadata.china$assembly_accession[i_idx], "_", metadata.china$assembly_accession[j_idx]),
  tip_i     = metadata.china$assembly_accession[i_idx],
  tip_j     = metadata.china$assembly_accession[j_idx],
  distance_km  = dist_km,
  stringsAsFactors = FALSE
)

metadata.pairwise <- metadata.pairwise %>% 
  dplyr::filter(distance_km != 0)

rm(n, pairs_mat, i_idx, j_idx, dist_km, dist_m)

# hist(metadata.pairwise$distance_km)



beast_tree <- treeio::read.beast('./data/beast_gubbins_contig/run_1/beast_GTRIG_clean_core_filtered.tree')

MAX_DATE <- 2023.58333333333
mrca.height <- beast_tree[, c("node", "height")]
mrca.height <- mrca.height %>% mutate(height = as.numeric(height))

# ---------------------------
# 生成两两样本对，并获取样本对的 MRCA 编号
# 
sample_ids <- intersect(metadata.china$assembly_accession, beast_tree@phylo$tip.label)
sample_pairs <- combn(sample_ids, 2)

# 获取 MRCA 编号
mrcas <- vapply(
  seq_len(ncol(sample_pairs)),
  function(i) getMRCA(beast_tree@phylo, sample_pairs[, i]),
  FUN.VALUE = integer(1)
)

mrca.nodes <- data.frame(
  sample_1 = sample_pairs[1, ],
  sample_2 = sample_pairs[2, ],
  mrca_node = mrcas
)

rm(sample_pairs, mrcas)

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
  ) %>% 
 mutate(year_diff = abs(collection_year_1 - collection_year_2))

mrca.nodes <- mrca.nodes %>% 
  inner_join(., dplyr::select(metadata.pairwise, -tip_pair), by = c("sample_1" = "tip_i", "sample_2" = "tip_j"))

mrca.nodes <- left_join(mrca.nodes, mrca.height, by = c('mrca_node' = 'node'))


# save.image('./r_image/distance_dynamic_resample_CI.RData')

################################################################################

load('./r_image/distance_dynamic_resample_CI.RData')


# 函数：给定数值向量，按分段点，生成有序分组
cut_into_group <- function(vec, breaks, unit = "", right_open = TRUE) {
  # vec: 需要分组的向量
  # breaks: 分段节点，比如 c(10, 50, 100, 500)
  # unit: 单位的间隔，比如 "y"（年）、"km"（公里）
  # right_open: 是否区间右开
  
  # 自动生成 labels
  labels <- c(
    paste0("<", breaks[1], unit),
    paste0(breaks[-length(breaks)], "-", breaks[-1], unit),
    paste0(">", breaks[length(breaks)], unit)
  )
  
  # 生成 cut breaks
  full_breaks <- c(-Inf, breaks, Inf)
  
  # 使用 cut 函数进行分组
  factor(
    cut(vec, breaks = full_breaks, labels = labels, right = right_open, ordered_result = TRUE),
    levels = labels
  )
}



mrca_cut_level <- c(10, 50, 100)
geo_cut_level  <- c(100, 200, 300, 400, 500, 1000, 2000, 3000)
ref_group <- "500-1000km"




mrca.nodes.l13 <- mrca.nodes %>%
  mutate(
    involves_lineage = (lineage_level_2_1 == "1.3" & geographic_province_1 == "Yunnan") | (lineage_level_2_2 == "1.3" & geographic_province_2 == "Yunnan")
  ) %>% 
  filter(involves_lineage)

ggplot(mrca.nodes.l13, aes(x = distance_km, y = height)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_point(color = 'navy', alpha = 0.2, size = 3) + 
  geom_smooth(fill = 'grey', linewidth = 1) + 
  theme_bw() + 
  annotation_logticks(sides = "bl")



mrca.nodes.l11 <- mrca.nodes %>%
  mutate(
    involves_lineage = (lineage_level_2_1 == "1.1" & geographic_province_1 == "Yunnan") | (lineage_level_2_2 == "1.1" & geographic_province_2 == "Yunnan")
  ) %>% 
  filter(involves_lineage) 

ggplot(mrca.nodes.l11, aes(x = distance_km, y = height)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_point(color = 'navy', alpha = 0.2, size = 3) + 
  geom_smooth(fill = 'grey', linewidth = 1) + 
  theme_bw() + 
  annotation_logticks(sides = "bl")



mrca.nodes.l12 <- mrca.nodes %>%
  mutate(
    involves_lineage = (lineage_level_2_1 == "1.2" & geographic_province_1 == "Yunnan") | (lineage_level_2_2 == "1.2" & geographic_province_2 == "Yunnan")
  ) %>% 
  filter(involves_lineage)

ggplot(mrca.nodes.l12, aes(x = distance_km, y = height)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_point(color = 'navy', alpha = 0.2, size = 3) + 
  geom_smooth(fill = 'grey', linewidth = 1) + 
  theme_bw() + 
  annotation_logticks(sides = "bl")



## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE: lineage_1.1_distance_RR_bootstrap_CI_centroid_points
mrca.nodes.fit <- mrca.nodes.l11 %>% 
  mutate(
    mrca_group = cut_into_group(height, breaks = mrca_cut_level, unit = "y"),
    geo_group  = cut_into_group(distance_km, breaks = geo_cut_level, unit = "km")
  ) %>%
  dplyr::filter(!is.na(geo_group) & !is.na(mrca_group) & year_diff <= 5)


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE: lineage_1.2_distance_RR_bootstrap_CI_centroid_points
mrca.nodes.fit <- mrca.nodes.l12 %>% 
  mutate(
    mrca_group = cut_into_group(height, breaks = mrca_cut_level, unit = "y"),
    geo_group  = cut_into_group(distance_km, breaks = geo_cut_level, unit = "km")
  ) %>%
  dplyr::filter(!is.na(geo_group) & !is.na(mrca_group) & year_diff <= 5)


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE: lineage_1.3_distance_RR_bootstrap_CI_centroid_points
mrca.nodes.fit <- mrca.nodes.l13 %>% 
  mutate(
    mrca_group = cut_into_group(height, breaks = mrca_cut_level, unit = "y"),
    geo_group  = cut_into_group(distance_km, breaks = geo_cut_level, unit = "km")
  ) %>%
  dplyr::filter(!is.na(geo_group) & !is.na(mrca_group) & year_diff <= 5)


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

table(mrca.nodes.fit$mrca_group, mrca.nodes.fit$geo_group)

# 设定需要保证组合完整的水平
all_mrca_groups <- levels(mrca.nodes.fit$mrca_group)
all_geo_groups  <- levels(mrca.nodes.fit$geo_group)

rr_summary <- mrca.nodes.fit %>%
  group_by(mrca_group, geo_group) %>%
  summarise(pair_count = n(), .groups = "drop") %>%
  # 补全所有组合，缺失的设为0
  complete(mrca_group = all_mrca_groups, geo_group = all_geo_groups, fill = list(pair_count = 0)) %>%
  # 把 pair_count = 0 的地方设为 1
  mutate(pair_count = ifelse(pair_count == 0, sample(seq(0.95, 1.05, 0.001), size = 1), pair_count)) %>%
  group_by(mrca_group) %>%
  mutate(
    mrca_group = factor(mrca_group, levels = all_mrca_groups), 
    geo_group = factor(geo_group, levels = all_geo_groups), 
    total = sum(pair_count),
    prob = pair_count / total,
    ref_prob = prob[geo_group == ref_group],
    RR = prob / ref_prob
  )

set.seed(1234)  # 确保可重复
n_bootstrap <- 1000

# 获取所有样本ID
all_samples <- unique(c(mrca.nodes.fit$sample_1, mrca.nodes.fit$sample_2))

# 保存每次bootstrap的RR结果
bootstrap_results <- vector("list", length = n_bootstrap)

for (b in seq_len(n_bootstrap)) {
  
  # (1) 有放回地bootstrap采样样本ID
  resampled_samples <- sample(all_samples, size = length(all_samples), replace = TRUE)
  resampled_samples <- unique(resampled_samples)  # 取unique以免组合太小
  
  # (2) 筛选pair：两个样本都在bootstrap集合中
  mrca.boot <- mrca.nodes.fit %>%
    filter(sample_1 %in% resampled_samples & sample_2 %in% resampled_samples)
  
  if (nrow(mrca.boot) == 0) next  # 避免空集合
  
  # (3) 计算每次bootstrap中的RR
  boot_rr <- mrca.boot %>%
    group_by(mrca_group, geo_group) %>%
    summarise(pair_count = n(), .groups = "drop") %>% 
    complete(mrca_group = all_mrca_groups, geo_group = all_geo_groups, fill = list(pair_count = 0)) %>% 
    mutate(pair_count = ifelse(pair_count == 0, sample(seq(0.95, 1.05, 0.001), size = 1), pair_count)) %>% 
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

# 合并所有bootstrap replicate结果
rr_boot_all <- bind_rows(bootstrap_results, .id = "replicate")  # .id 是bootstrap编号

# 计算每个(mrca_group, geo_group)对应的2.5%、97.5%分位数
rr_ci <- rr_boot_all %>%
  group_by(mrca_group, geo_group) %>%
  summarise(
    Meddian  = quantile(RR, probs = 0.5, na.rm = TRUE),
    RR_lower = quantile(RR, probs = 0.025, na.rm = TRUE),
    RR_upper = quantile(RR, probs = 0.975, na.rm = TRUE), 
    RR_lower_inter = Meddian - RR_lower, 
    RR_upper_inter = RR_upper - Meddian, 
    .groups = "drop"
  ) %>% 
  dplyr::select(-c(RR_lower, RR_upper))

# 合并回原rr_summary
rr_summary <- rr_summary %>%
  left_join(rr_ci, by = c("mrca_group", "geo_group")) %>% 
  mutate(
    RR_lower = RR - RR_lower_inter, 
    RR_upper = RR + RR_upper_inter, 
    mrca_group = factor(mrca_group, levels = all_mrca_groups), 
    geo_group = factor(geo_group, levels = all_geo_groups), 
  ) %>% 
  dplyr::select(-c(RR_lower_inter, RR_upper_inter))

# write.csv(rr_summary, "./result/geospatial_dynamic/lineage_1.1_distance_RR_bootstrap_CI_centroid_points_summzry.csv", row.names = FALSE)
# write.csv(rr_summary, "./result/geospatial_dynamic/lineage_1.2_distance_RR_bootstrap_CI_centroid_points_summzry.csv", row.names = FALSE)
# write.csv(rr_summary, "./result/geospatial_dynamic/lineage_1.3_distance_RR_bootstrap_CI_centroid_points_summzry.csv", row.names = FALSE)


ggplot(rr_summary, aes(x = RR, y = geo_group)) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +  # 参考线
  geom_errorbarh(aes(xmin = RR_lower, xmax = RR_upper), height = 0.1, color = "black", linewidth = 0.5) +  # 置信区间
  geom_point(size = 2.5, color = "black") +  # RR 点
  scale_x_log10() + 
  theme_classic(base_size = 14) + 
  coord_flip() + 
  facet_wrap(~ mrca_group, ncol = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



################################################################################

rm(list = ls())
gc()

rr_l11 <-  read.csv("./result/geospatial_dynamic/lineage_1.1_distance_RR_bootstrap_CI_centroid_points_summzry.csv")
rr_l12 <-  read.csv("./result/geospatial_dynamic/lineage_1.2_distance_RR_bootstrap_CI_centroid_points_summzry.csv")
rr_l13 <-  read.csv("./result/geospatial_dynamic/lineage_1.3_distance_RR_bootstrap_CI_centroid_points_summzry.csv")


rr_main <- rbind(
  rr_l11 %>% mutate(Lineage = "1.1"), 
  rr_l12 %>% mutate(Lineage = "1.2"), 
  rr_l13 %>% mutate(Lineage = "1.3")
) %>% 
  dplyr::filter(!(geo_group != "500-1000km" & RR == 1) & RR_lower > 0) %>% 
  mutate(
    mrca_group = factor(mrca_group, levels = c("<10y", "10-50y", "50-100y", ">100y")), 
    geo_group  = factor(geo_group,  levels = c("<100km", "100-200km", "200-300km", "300-400km", "400-500km", "500-1000km", "1000-2000km", "2000-3000km", ">3000km"))
  )

ggplot(rr_main, aes(x = geo_group, y = RR)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +  # 参考线
  geom_errorbar(aes(ymin = RR_lower, ymax = RR_upper, color = Lineage), width = 0.2, linewidth = 0.5) +  # 置信区间 , color = "black"
  # geom_point(size = 2.5, color = "black") +  # RR 点
  geom_point(data = dplyr::filter(rr_main, geo_group != "500-1000km"), aes(color = Lineage), size = 3, shape = 17) +
  geom_point(data = dplyr::filter(rr_main, geo_group == "500-1000km"), aes(color = Lineage), size = 3, shape = 0, stroke = 1) +
  scale_y_log10() + 
  theme_classic(base_size = 14) + 
  facet_wrap(mrca_group ~ Lineage, ncol = 3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(
    x = "Geographic distance (km)", 
    y = "Relative Risk"
  ) + 
  scale_color_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                                "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                                "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                     na.value = "transparent")

# ggsave('./plot/geospatial_dynamic/all_lineage_and_mrca_group_distance_dynamic_resample_CI_centroid_points_plot.pdf', width = 9, height = 12)

