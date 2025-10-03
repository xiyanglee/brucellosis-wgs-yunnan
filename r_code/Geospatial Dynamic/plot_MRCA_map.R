library(cowplot)
library(tidyverse)
library(sf)

load('./r_image/distance_dynamic_random_points.RData')

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')



metadata.china <- metadata.china %>% 
  left_join(., metadata.phylo[, c("assembly_accession", "lineage_level_2")], by = c("assembly_accession"))
metadata.china <- metadata.china[!is.na(metadata.china$lineage_level_2), ]

metadata.china %>% 
  dplyr::filter(lineage_level_2 == "1.3" & lng > 90 & lng < 100 & lat > 30 & lat < 40)

metadata.china[which(metadata.china$assembly_accession == "GCA_032252615.1"), "lng"] <- 110.8781 # 103.8781
metadata.china[which(metadata.china$assembly_accession == "GCA_032252615.1"), "lat"] <- 41.90776 # 40.90776
metadata.china[which(metadata.china$assembly_accession == "GCA_016806105.1"), "lng"] <- 96.7922 # 98.89221
metadata.china[which(metadata.china$assembly_accession == "GCA_016806105.1"), "lat"] <- 36.44443 # 36.34443

mrca.plot <- mrca.nodes %>% 
  # dplyr::filter(height < 100) %>% 
  mutate(
    involves_lineage = (geographic_province_1 == "Yunnan") | (geographic_province_2 == "Yunnan")
  ) %>% 
  filter(involves_lineage)
  
mrca.plot <- mrca.plot %>% 
  left_join(
    ., 
    metadata.china %>% 
      dplyr::select(assembly_accession, lng, lat) %>%
      rename_with(~paste0(., "_1"), .cols = c(lng, lat)), 
    by = c("sample_1" = "assembly_accession")
  ) %>% 
  left_join(
    ., 
    metadata.china %>% 
      dplyr::select(assembly_accession, lng, lat) %>%
      rename_with(~paste0(., "_2"), .cols = c(lng, lat)), 
    by = c("sample_2" = "assembly_accession")
  )



# 加载省级地图
geo_sf.province <- read_sf('../../metadata/中国_省.geojson')

# 方位投影字符串（Azimuthal Equidistant Projection）
proj_aeqd <- "+proj=aeqd +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# 准备 mrca.plot 中的样本对连线数据，转为 sf 的 LINESTRING
lines_sf <- mrca.plot %>%
  rowwise() %>%
  mutate(geometry = st_sfc(
    st_linestring(matrix(c(lng_1, lat_1, lng_2, lat_2), ncol = 2, byrow = TRUE)),
    crs = 4326
  )) %>%
  st_as_sf()

# 准备 metadata.china 的样本点数据
points_sf <- metadata.china %>%
  st_as_sf(coords = c("lng", "lat"), crs = 4326)

# 将所有图层都统一投影到 aeqd 投影
geo_sf.province_proj <- st_transform(geo_sf.province, crs = proj_aeqd)
lines_sf_proj <- st_transform(lines_sf, crs = proj_aeqd)
points_sf_proj <- st_transform(points_sf, crs = proj_aeqd)
# 中国外边界轮廓
geo_sf.china_outline <- geo_sf.province_proj %>% st_union()
## 确定三类地区
## 农业农村部关于印发《国家动物疫病监测与流行病学调查计划（2021—2025年）》的通知
## https://www.moa.gov.cn/nybgb/2021/202105/202110/t20211021_6380181.htm
## 一类地区: 人间报告发病率超过十万分之一或畜间疫情未控制县数占总县数30%以上的省份，包括北京、天津、河北、山西、内蒙古、辽宁、吉林、黑龙江、山东、河南、陕西、甘肃、青海、宁夏、新疆等15个省份和新疆生产建设兵团
## 二类地区: 本地有新发人间病例发生且报告发病率低于或等于十万分之一、畜间疫情未控制县数占总县数30%以下的省份，包括上海、江苏、浙江、安徽、福建、江西、湖北、湖南、广东、广西、重庆、四川、贵州、云南、西藏等15个省份
## 三类地区: 无本地新发人间病例和畜间疫情的省份，目前有海南省
geo_sf.province_proj <- geo_sf.province_proj %>%
  mutate(region_class = case_when(
    name %in% c('北京市', '天津市', '河北省', '山西省', '内蒙古自治区', '辽宁省', '吉林省', '黑龙江省', '山东省', '河南省', '陕西省', '甘肃省', '青海省', '宁夏回族自治区', '新疆维吾尔自治区') ~ 'Class I',
    name %in% c('上海市', '江苏省', '浙江省', '安徽省', '福建省', '江西省', '湖北省', '湖南省', '广东省', '广西壮族自治区', '重庆市', '四川省', '贵州省', '云南省', '西藏自治区') ~ 'Class II',
    name %in% c('海南省') ~ 'Class III',
    .default = 'Others'
  ))

main_map <- ggplot() + 
  geom_sf(data = geo_sf.province_proj, fill = "white", color = "black", linewidth = 0.3) + 
  geom_sf(data = dplyr::filter(lines_sf_proj, height >= 100 | year_diff > 5), color = "#f1f2ed", size = 1, alpha = 0.8) + 
  geom_sf(data = dplyr::filter(lines_sf_proj, height <  100 & height >= 10 & year_diff <= 5), aes(color = height), size = 1, alpha = 0.8) +
  geom_sf(data = dplyr::filter(lines_sf_proj, height <= 15  & year_diff <= 5), aes(color = height), size = 1, alpha = 0.8) +
  scale_color_gradientn(
    name = "Height",
    colors = c("#007991", "#007991", "#e0f2b6", "#eff2e6"),
    values = scales::rescale(c(0, 5, 25, 50, 100), to = c(0, 1)),
    limits = c(0, 100),
    na.value = "grey50"
  ) + 
  geom_sf(data = points_sf_proj, aes(fill = lineage_level_2), color = "black", size = 3, shape = 21) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                                "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                                "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                     na.value = "transparent") + 
  theme_void() +
  theme(
    legend.position = "right",
    # panel.background = element_rect(fill = "white"),
    # panel.grid = element_line(color = "grey90")
  ) + 
  coord_sf(
    crs = proj_aeqd,   # 投影
    xlim = c(70, 140), # 经纬度
    ylim = c(18, 55),
    default_crs = sf::st_crs(4326), # 指示 xlim / ylim 是 EPSG:4326 坐标
    expand = FALSE
  )

south_china_sea_map <- ggplot() + 
  geom_sf(data = geo_sf.province_proj, fill = "white", color = "black", linewidth = 0.3) +
  theme_void() + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  ) + 
  coord_sf(
    crs = proj_aeqd, 
    xlim = c(106.1, 121), 
    ylim = c(0, 23.5),
    default_crs = sf::st_crs(4326), 
    expand = FALSE
  )

ggdraw() +
  draw_plot(main_map) +
  draw_plot(south_china_sea_map, x = 0.675, y = 0.1, width = 0.25, height = 0.275)

# ggsave('./plot/geospatial_dynamic/MRCA_map.pdf', width = 9, height = 7)



################################################################################
## Lineage 1.1

points_sf_proj_target <- points_sf_proj %>% 
  dplyr::filter(lineage_level_2 == "1.3" & str_starts(assembly_accession, "20"))

points_sf_proj_others <- points_sf_proj %>% 
  dplyr::filter( ! assembly_accession %in% points_sf_proj_target$assembly_accession)

# 创建矩形框（原始 CRS = WGS84）
box_proj <- st_as_sfc(st_bbox(c(xmin = 98.5, xmax = 105.1, ymin = 22.2, ymax = 28.1), crs = 4326))
# 投影为地图当前的 CRS（例如 aeqd）
box_proj <- st_transform(box_proj, crs = st_crs(geo_sf.province_proj))

# min(lines_sf_proj$height)
lines_sf_proj_subset <- lines_sf_proj %>% 
  mutate(
    involves_lineage = (lineage_level_2_1 == "1.1" & geographic_province_1 == "Yunnan") | (lineage_level_2_2 == "1.1" & geographic_province_2 == "Yunnan")
  ) %>% 
  filter(involves_lineage) %>% 
  # ungroup() + 
  mutate(
    height = height - 5, # min height: 2.972263
    mrca_group = case_when(
      height < 20   ~ "0-10 y",
      height >= 20  & height < 50  ~ "10-50 y",
      height >= 50  & height < 100 ~ "50-100 y",
      height >= 100                ~ ">100 y",
      TRUE                         ~ NA_character_
    )
  )

main_map <- ggplot() + 
  # 地图
  
  # geom_sf(data = dplyr::filter(geo_sf.province_proj, region_class == "Class I"), fill = "#f5efe9", linewidth = 0.3, color = "white", alpha = 1) +
  # geom_sf(data = dplyr::filter(geo_sf.province_proj, region_class == "Class II"), fill = "#f0f3f5", linewidth = 0.3, color = "white", alpha = 1) +
  # geom_sf(data = geo_sf.china_outline, fill = NA, linewidth = 0.3, color = "black") +
  
  geom_sf(data = geo_sf.province_proj, fill = "white", linewidth = 0.3, color = "grey90") +
  geom_sf(data = geo_sf.china_outline, fill = NA, linewidth = 0.3, color = "black") +
  
  # 样本对连线，分组以控制涂层顺序
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height >= 100 | year_diff > 5), aes(color = mrca_group), size = 1, alpha = 0.5) + 
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height <  100 & height >= 10 & year_diff <= 5), aes(color = mrca_group), size = 1, alpha = 0.8) +
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height <= 15  & year_diff <= 5), aes(color = mrca_group), size = 1.5, alpha = 1) +
  # scale_color_gradientn(
  #   name = "Height",
  #   colors = c("#007991", "#007991", "#e0f2b6", "#eff2e6"),
  #   values = scales::rescale(c(0, 10, 20, 50, 100), to = c(0, 1)),
  #   limits = c(0, 100),
  #   na.value = "grey50"
  # ) +
  scale_color_manual(name = "tMRCA", values = c("0-10 y" = "#007991", "10-50 y" = "#99c4cc", "50-100 y" = "#e0f2b6", ">100 y" = "#f0f2eb")) +
  
  # 样本点
  geom_sf(data = points_sf_proj, fill = NA, color = "black", size = 2, shape = 21, alpha = 1) + 
  geom_sf(data = points_sf_proj_others, aes(fill = lineage_level_2), color = "black", size = 2, shape = 21, alpha = 0.8) + 
  geom_sf(data = points_sf_proj_target, aes(fill = lineage_level_2), color = "black", size = 2, shape = 21, alpha = 0.8) + 
  # geom_sf(data = points_sf_proj, aes(fill = lineage_level_2), color = "black", size = 3, shape = 21, alpha = 0.8) + 
  scale_fill_manual(name = "Lineage", 
                    values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                               "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                               "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                    na.value = "transparent") + 
  
  
  theme_void() +
  theme(legend.position = "right") + 
  coord_sf(
    crs = proj_aeqd,   # 投影
    # xlim = c(70, 140), # 经纬度
    # ylim = c(18, 55),
    # default_crs = sf::st_crs(4326), # 指示 xlim / ylim 是 EPSG:4326 坐标
    expand = FALSE
  ) + 
  geom_sf(data = box_proj, fill = NA, color = "black", linewidth = 0.5, linetype = "dashed")

yunnan_map <- main_map + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), 
    legend.position = "none"
  ) + 
  coord_sf(
    crs = proj_aeqd,   # 投影
    xlim = c(98.5, 105.1), # 经纬度
    ylim = c(22.2, 28.1),
    default_crs = sf::st_crs(4326), # 指示 xlim / ylim 是 EPSG:4326 坐标
    expand = FALSE
  )

ggdraw() + 
  draw_plot(main_map) +
  draw_plot(yunnan_map, x = 0.1, y = 0.075, width = 0.3, height = 0.3)

# ggsave('./plot/geospatial_dynamic/l11_MRCA_map.pdf', width = 5.5, height = 5)

################################################################################
## Lineage 1.2

points_sf_proj_target <- points_sf_proj %>% 
  dplyr::filter(lineage_level_2 == "1.2" & str_starts(assembly_accession, "20"))

points_sf_proj_others <- points_sf_proj %>% 
  dplyr::filter( ! assembly_accession %in% points_sf_proj_target$assembly_accession)

# 创建矩形框（原始 CRS = WGS84）
box_proj <- st_as_sfc(st_bbox(c(xmin = 98.5, xmax = 105.1, ymin = 22.2, ymax = 28.1), crs = 4326))
# 投影为地图当前的 CRS（例如 aeqd）
box_proj <- st_transform(box_proj, crs = st_crs(geo_sf.province_proj))

# min(lines_sf_proj$height)
lines_sf_proj_subset <- lines_sf_proj %>% 
  mutate(
    involves_lineage = (lineage_level_2_1 == "1.2" & geographic_province_1 == "Yunnan") | (lineage_level_2_2 == "1.2" & geographic_province_2 == "Yunnan")
  ) %>% 
  filter(involves_lineage) %>% 
  # ungroup() + 
  mutate(
    height = height - 5, # min height: 2.972263
    mrca_group = case_when(
      height < 20   ~ "0-10 y",
      height >= 20  & height < 50  ~ "10-50 y",
      height >= 50  & height < 100 ~ "50-100 y",
      height >= 100                ~ ">100 y",
      TRUE                         ~ NA_character_
    )
  )

main_map <- ggplot() + 
  # 地图
  
  # geom_sf(data = dplyr::filter(geo_sf.province_proj, region_class == "Class I"), fill = "#f5efe9", linewidth = 0.3, color = "white", alpha = 1) +
  # geom_sf(data = dplyr::filter(geo_sf.province_proj, region_class == "Class II"), fill = "#f0f3f5", linewidth = 0.3, color = "white", alpha = 1) +
  # geom_sf(data = geo_sf.china_outline, fill = NA, linewidth = 0.3, color = "black") +
  
  geom_sf(data = geo_sf.province_proj, fill = "white", linewidth = 0.3, color = "grey90") +
  geom_sf(data = geo_sf.china_outline, fill = NA, linewidth = 0.3, color = "black") +
  
  # 样本对连线，分组以控制涂层顺序
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height >= 100 | year_diff > 5), aes(color = mrca_group), size = 1, alpha = 0.5) + 
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height <  100 & height >= 10 & year_diff <= 5), aes(color = mrca_group), size = 1, alpha = 0.8) +
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height <= 15  & year_diff <= 5), aes(color = mrca_group), size = 1.5, alpha = 1) +
  scale_color_manual(name = "tMRCA", values = c("0-10 y" = "#007991", "10-50 y" = "#99c4cc", "50-100 y" = "#e0f2b6", ">100 y" = "#f0f2eb")) +
  
  # 样本点
  geom_sf(data = points_sf_proj, fill = NA, color = "black", size = 2, shape = 21, alpha = 1) + 
  geom_sf(data = points_sf_proj_others, aes(fill = lineage_level_2), color = "black", size = 2, shape = 21, alpha = 0.8) + 
  geom_sf(data = points_sf_proj_target, aes(fill = lineage_level_2), color = "black", size = 2, shape = 21, alpha = 0.8) + 
  scale_fill_manual(name = "Lineage", 
                    values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                               "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                               "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                    na.value = "transparent") + 
  
  
  theme_void() +
  theme(legend.position = "right") + 
  coord_sf(
    crs = proj_aeqd,   # 投影
    expand = FALSE
  ) + 
  geom_sf(data = box_proj, fill = NA, color = "black", linewidth = 0.5, linetype = "dashed")

yunnan_map <- main_map + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), 
    legend.position = "none"
  ) + 
  coord_sf(
    crs = proj_aeqd,   # 投影
    xlim = c(98.5, 105.1), # 经纬度
    ylim = c(22.2, 28.1),
    default_crs = sf::st_crs(4326), # 指示 xlim / ylim 是 EPSG:4326 坐标
    expand = FALSE
  )

ggdraw() + 
  draw_plot(main_map) +
  draw_plot(yunnan_map, x = 0.1, y = 0.075, width = 0.3, height = 0.3)

# ggsave('./plot/geospatial_dynamic/l12_MRCA_map.pdf', width = 5.5, height = 5)

################################################################################
## Lineage 1.3

points_sf_proj_target <- points_sf_proj %>% 
  dplyr::filter(lineage_level_2 == "1.3" & str_starts(assembly_accession, "20"))

points_sf_proj_others <- points_sf_proj %>% 
  dplyr::filter( ! assembly_accession %in% points_sf_proj_target$assembly_accession)

# 创建矩形框（原始 CRS = WGS84）
box_proj <- st_as_sfc(st_bbox(c(xmin = 98.5, xmax = 105.1, ymin = 22.2, ymax = 28.1), crs = 4326))
# 投影为地图当前的 CRS（例如 aeqd）
box_proj <- st_transform(box_proj, crs = st_crs(geo_sf.province_proj))

# min(lines_sf_proj$height)
lines_sf_proj_subset <- lines_sf_proj %>% 
  mutate(
    involves_lineage = (lineage_level_2_1 == "1.3" & geographic_province_1 == "Yunnan") | (lineage_level_2_2 == "1.3" & geographic_province_2 == "Yunnan")
  ) %>% 
  filter(involves_lineage) %>% 
  # ungroup() + 
  mutate(
    height = height - 5, # min height: 2.972263
    mrca_group = case_when(
      height < 20   ~ "0-10 y",
      height >= 20  & height < 50  ~ "10-50 y",
      height >= 50  & height < 100 ~ "50-100 y",
      height >= 100                ~ ">100 y",
      TRUE                         ~ NA_character_
    )
  )

main_map <- ggplot() + 
  # 地图
  
  # geom_sf(data = dplyr::filter(geo_sf.province_proj, region_class == "Class I"), fill = "#f5efe9", linewidth = 0.3, color = "white", alpha = 1) +
  # geom_sf(data = dplyr::filter(geo_sf.province_proj, region_class == "Class II"), fill = "#f0f3f5", linewidth = 0.3, color = "white", alpha = 1) +
  # geom_sf(data = geo_sf.china_outline, fill = NA, linewidth = 0.3, color = "black") +
  
  geom_sf(data = geo_sf.province_proj, fill = "white", linewidth = 0.3, color = "grey90") +
  geom_sf(data = geo_sf.china_outline, fill = NA, linewidth = 0.3, color = "black") +
  
  # 样本对连线，分组以控制涂层顺序
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height >= 100 | year_diff > 5), aes(color = mrca_group), size = 1, alpha = 0.5) + 
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height <  100 & height >= 10 & year_diff <= 5), aes(color = mrca_group), size = 1, alpha = 0.8) +
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height <= 15  & year_diff <= 5), aes(color = mrca_group), size = 1.5, alpha = 1) +
  # scale_color_gradientn(
  #   name = "Height",
  #   colors = c("#007991", "#007991", "#e0f2b6", "#eff2e6"),
  #   values = scales::rescale(c(0, 10, 20, 50, 100), to = c(0, 1)),
  #   limits = c(0, 100),
  #   na.value = "grey50"
  # ) +
  scale_color_manual(values = c("0-10 y" = "#007991", "10-50 y" = "#99c4cc", "50-100 y" = "#e0f2b6", ">100 y" = "#f0f2eb")) +
  
  # 样本点
  geom_sf(data = points_sf_proj, fill = NA, color = "black", size = 2, shape = 21, alpha = 1) + 
  geom_sf(data = points_sf_proj_others, aes(fill = lineage_level_2), color = "black", size = 2, shape = 21, alpha = 0.8) + 
  ##############################################################################
  ## 再次绘制 tMRCA 较小的连接线 
  geom_sf(data = dplyr::filter(lines_sf_proj_subset, height <= 15  & year_diff <= 5), aes(color = mrca_group), size = 1.5, alpha = 1) + 
  geom_sf(data = points_sf_proj_target, aes(fill = lineage_level_2), color = "black", size = 2, shape = 21, alpha = 0.8) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                               "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                               "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                    na.value = "transparent") + 
  
  
  theme_void() +
  theme(legend.position = "right") + 
  coord_sf(
    crs = proj_aeqd,   # 投影
    # xlim = c(70, 140), # 经纬度
    # ylim = c(18, 55),
    # default_crs = sf::st_crs(4326), # 指示 xlim / ylim 是 EPSG:4326 坐标
    expand = FALSE
  ) + 
  geom_sf(data = box_proj, fill = NA, color = "black", linewidth = 0.5, linetype = "dashed")

yunnan_map <- main_map + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), 
    legend.position = "none"
  ) + 
  coord_sf(
    crs = proj_aeqd,   # 投影
    xlim = c(98.5, 105.1), # 经纬度
    ylim = c(22.2, 28.1),
    default_crs = sf::st_crs(4326), # 指示 xlim / ylim 是 EPSG:4326 坐标
    expand = FALSE
  )

ggdraw() + 
  draw_plot(main_map) +
  draw_plot(yunnan_map, x = 0.1, y = 0.075, width = 0.3, height = 0.3)

# ggsave('./plot/geospatial_dynamic/l13_MRCA_map.pdf', width = 5.5, height = 5)

