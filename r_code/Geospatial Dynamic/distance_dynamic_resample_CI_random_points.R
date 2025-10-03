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

## 海南省范围包括南海海域
## 由于之后需要对各省份样本进行随机定点，包括海域大概率会抽到南海某点
## 计算面积最大多边形 (海南主岛) 并替换
hainan_index <- which(geo_sf.province$name == "海南省")
hainan_parts <- geo_sf.province[hainan_index, ]
# 只提取 geometry，变成 sfc（Simple Feature Collection）
hainan_parts <- st_geometry(hainan_parts)
# 将 MULTIPOLYGON 拆分成 POLYGON
hainan_parts <- st_cast(hainan_parts, "POLYGON")
# 变回 sf 对象
hainan_parts <- st_sf(geometry = hainan_parts)
# 计算每个 polygon 的面积
hainan_parts <- hainan_parts %>% mutate(area = st_area(geometry))
# 选择面积最大的 polygon （即陆地）
hainan_parts <- hainan_parts %>%
  filter(area == max(area)) %>%
  select(-area)
# 4. 替换 geo_sf.province 中海南省的 geometry
geo_sf.province$geometry[hainan_index] <- hainan_parts$geometry
rm(hainan_index, hainan_parts)

geo_sf <- bind_rows(geo_sf.province, geo_sf.city, geo_sf.county) %>% distinct()
rm(geo_sf.province, geo_sf.city, geo_sf.county)

head(sort(table(geo_sf$name), decreasing = TRUE), 32)
# 境界线 市中区 鼓楼区   城区 新华区   郊区 铁西区   东区 南山区 南沙区 向阳区 
#     15      4      4      3      3      3      3      2      2      2      2 
# 和平区 城中区 城关区 宝山区 新城区 普陀区 朝阳区 桥西区 永定区 江北区 河东区 
#      2      2      2      2      2      2      2      2      2      2      2 
# 海州区 白云区 西安区 西湖区 通州区 铁东区 长安区 青山区 龙华区 丁青县 
#      2      2      2      2      2      2      2      2      2      1 

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


# 使用 apply 函数进行随机抽样
random_points <- apply(metadata.china, 1, function(data) {
  st_sample(data$geometry, size = 1)  # 从每个几何区域抽取一个点
})

# 提取经度和纬度
lng <- sapply(random_points, function(point) st_coordinates(point)[1])  # 提取经度
lat <- sapply(random_points, function(point) st_coordinates(point)[2])  # 提取纬度

# 将经度和纬度添加到原始数据中
metadata.china$lng <- lng
metadata.china$lat <- lat

metadata.china <- metadata.china %>% 
  select(-c(gb, geometry)) %>% 
  st_drop_geometry()

rm(lng, lat, random_points, geo_sf)

# plot(metadata.china$lng, metadata.china$lat)


# 取两个表格列的交集并合并
common_cols    <- intersect(colnames(metadata.china), colnames(metadata.sample))
metadata.china <- bind_rows(
  metadata.china[, common_cols],
  metadata.sample[, common_cols]
)
rm(common_cols)

# plot(metadata.china$lng, metadata.china$lat)

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

# save.image('./r_image/distance_dynamic_random_points.RData')

################################################################################

load('./r_image/distance_dynamic_random_points.RData')

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

head(sort(table(mrca.nodes.l12$height), decreasing = TRUE))

ggplot(mrca.nodes.l12, aes(x = distance_km, y = height)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_point(color = 'navy', alpha = 0.2, size = 3) + 
  geom_smooth(fill = 'grey', linewidth = 1) + 
  theme_bw() + 
  annotation_logticks(sides = "bl")



## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE: lineage_1.1_distance_RR_bootstrap_CI_random_points
mrca.nodes.fit <- mrca.nodes.l11 %>% 
  mutate(
    mrca_group = cut_into_group(height, breaks = mrca_cut_level, unit = "y"),
    geo_group  = cut_into_group(distance_km, breaks = geo_cut_level, unit = "km")
  ) %>%
  dplyr::filter(!is.na(geo_group) & !is.na(mrca_group) & year_diff <= 5)


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE: lineage_1.3_distance_RR_bootstrap_CI_random_points
mrca.nodes.fit <- mrca.nodes.l13 %>% 
  mutate(
    mrca_group = cut_into_group(height, breaks = mrca_cut_level, unit = "y"),
    geo_group  = cut_into_group(distance_km, breaks = geo_cut_level, unit = "km")
  ) %>%
  dplyr::filter(!is.na(geo_group) & !is.na(mrca_group) & year_diff <= 5)


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE: lineage_1.2_distance_RR_bootstrap_CI_random_points
mrca.nodes.fit <- mrca.nodes.l12 %>% 
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

# write.csv(rr_summary, "./result/geospatial_dynamic/lineage_1.1_distance_RR_bootstrap_CI_random_points_summzry.csv", row.names = FALSE)
# write.csv(rr_summary, "./result/geospatial_dynamic/lineage_1.2_distance_RR_bootstrap_CI_random_points_summzry.csv", row.names = FALSE)
# write.csv(rr_summary, "./result/geospatial_dynamic/lineage_1.3_distance_RR_bootstrap_CI_random_points_summzry.csv", row.names = FALSE)


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

rr_l11 <-  read.csv("./result/geospatial_dynamic/lineage_1.1_distance_RR_bootstrap_CI_random_points_summzry.csv")
rr_l12 <-  read.csv("./result/geospatial_dynamic/lineage_1.2_distance_RR_bootstrap_CI_random_points_summzry.csv")
rr_l13 <-  read.csv("./result/geospatial_dynamic/lineage_1.3_distance_RR_bootstrap_CI_random_points_summzry.csv")

rr_main <- rbind(
  ## Lineage 1.2 单独画图
  # rr_l12 %>% dplyr::filter(mrca_group == "10-50y") %>% mutate(Lineage = "1.2")
  ## Lineage 1.1 + 1.3 合并画图
  rr_l11 %>% dplyr::filter(mrca_group == "10-50y") %>% mutate(Lineage = "1.1"),
  rr_l13 %>% dplyr::filter(mrca_group == "10-50y") %>% mutate(Lineage = "1.3")
) %>% 
  mutate(
    geo_distance = case_when(
      geo_group == "<100km" ~ 50, 
      geo_group == "100-200km" ~ 150, 
      geo_group == "200-300km" ~ 250, 
      geo_group == "300-400km" ~ 350, 
      geo_group == "400-500km" ~ 450, 
      geo_group == "500-1000km" ~ 750, 
      geo_group == "1000-2000km" ~ 1500, 
      geo_group == "2000-3000km" ~ 2500, 
      geo_group == ">3000km" ~ 3500
    )
  ) %>% 
  dplyr::filter(!(geo_group != "500-1000km" & RR == 1))

# ggplot(rr_main, aes(x = RR, y = geo_distance)) + 
#   geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +  # 参考线
#   geom_errorbarh(aes(xmin = RR_lower, xmax = RR_upper, color = Lineage), height = 0.1, linewidth = 0.5) +  # 置信区间
#   geom_point(aes(color = Lineage), size = 2.5) +  # RR 点
#   scale_x_log10() + 
#   scale_y_sqrt() + 
#   theme_classic(base_size = 14) + 
#   coord_flip() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(rr_main, aes(y = RR, x = geo_distance)) + 
#   # 绘制背景色
#   annotate("rect", xmin = 50, xmax = 750, ymin = 1, ymax = 100, fill = "#e6e4e3") + 
#   annotate("rect", xmin = 750, xmax = 3000, ymin = 1, ymax = 100, fill = "#edf5e6") + 
#   geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +  # 参考线
#   # geom_ribbon(aes(ymin = RR_lower, ymax = RR_upper, fill = Lineage), alpha = 0.3) +  # 使用 geom_ribbon 代替 geom_errorbarh
#   # geom_errorbar(aes(ymin = RR_lower, ymax = RR_upper, color = Lineage), width = 0, linewidth = 0.5) +
#   geom_line(aes(y = RR_upper, color = Lineage), linewidth = 0.5, linetype = "dashed") +
#   geom_line(aes(y = RR_lower, color = Lineage), linewidth = 0.5, linetype = "dashed") +
#   geom_line(aes(color = Lineage), linewidth = 1) +  # 使用 geom_line 代替 geom_point
#   geom_point(aes(color = Lineage), size = 2.5) +  # RR 点
#   # scale_x_sqrt(breaks = c(100, 200, 300, 400, 500, 1000, 2000, 3000)) +  # 设置 x 轴标签
#   scale_x_log10(breaks = c(100, 200, 300, 500, 1000, 2000, 3000)) + 
#   scale_y_log10() +  # 使用 log10 变换 y 轴
#   annotation_logticks(sides = "b") +
#   theme_classic(base_size = 14) + 
#   facet_wrap(~ Lineage, ncol = 2) + 
#   theme(strip.text = element_blank()) + # 不显示分面标题 
#   scale_color_manual(values = c("1.1" = "#7EBB8F", "1.3" = "#886A9E"), na.value = "transparent") + 
#   scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.3" = "#886A9E"), na.value = "transparent") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 设置 x 轴文本旋转

ggplot(rr_main, aes(y = RR, x = geo_distance)) + 
  # 背景色
  annotate("rect", xmin = 45, xmax = 500, ymin = 0.5, ymax = Inf, fill = "#f2f5f0") + 
  annotate("rect", xmin = 500, xmax = 1000, ymin = 0.5, ymax = Inf, fill = "grey98") + 
  annotate("rect", xmin = 1000, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "#f5f2f0") + 
  # geom_vline(xintercept = 750, linetype = "dashed", color = "grey40", linewidth = 0.8) + 
  
  # 主图
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) + 
  geom_ribbon(aes(ymin = RR_lower, ymax = RR_upper, fill = Lineage), alpha = 0.1) + 
  geom_line(aes(y = RR_upper, color = Lineage), linewidth = 0.5, linetype = "dashed") +
  geom_line(aes(y = RR_lower, color = Lineage), linewidth = 0.5, linetype = "dashed") +
  geom_line(aes(color = Lineage), linewidth = 1) +
  geom_point(data = dplyr::filter(rr_main, geo_distance != 750), aes(color = Lineage), size = 3, shape = 17) + 
  geom_point(data = dplyr::filter(rr_main, geo_distance == 750), aes(color = Lineage), size = 3, shape = 0, stroke = 1) + 
  
  # 坐标轴
  scale_x_log10(
    breaks = c(100, 200, 300, 500, 1000, 2000, 3000), 
    limits = c(45, 3800), 
    expand = c(0, 0)
  ) + 
  scale_y_log10(
    expand = c(0, 0)
  ) + 
  annotation_logticks(sides = "b") +
  
  # 样式
  theme_classic(base_size = 14) + 
  facet_wrap(~ Lineage, ncol = 2) + 
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + 
  labs(
    x = "Geographic distance (km)",
    y = "Relative Risk"
  ) + 
  
  # 配色
  scale_color_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E"), na.value = "transparent") + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E"), na.value = "transparent")

# ggsave('./plot/geospatial_dynamic/l11_l13_distance_dynamic_resample_CI_random_points_plot.pdf', width = 9.5, height = 3.5)
# ggsave('./plot/geospatial_dynamic/l12_distance_dynamic_resample_CI_random_points_plot.pdf', width = 5, height = 3.5)

write.csv(rr_main, "./r_plot_data/figure_5c_rr_main.csv", row.names = FALSE)

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

# ggsave('./plot/geospatial_dynamic/all_lineage_and_mrca_group_distance_dynamic_resample_CI_random_points_plot.pdf', width = 9, height = 12)

