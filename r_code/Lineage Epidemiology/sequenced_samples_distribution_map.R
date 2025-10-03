# install.packages('spatstat', lib = '/public/r_share_library/4.1')
# install.packages('rnaturalearth', lib = '/public/r_share_library/4.1')
# install.packages('rnaturalearthdata', lib = '/public/r_share_library/4.1')

library(terra)      # 栅格数据处理和制图
library(tidyterra)  # 栅格数据处理和制图
library(sf)         # 处理矢量数据
library(ggplot2)    # 创建地图
library(cowplot)    # 结合ggplot2对象
library(ggspatial)  # 添加指北针
library(elevatr)    # 用于下载 DEM 高程数据
library(spatstat)   # 计算空间密度
library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v2.RData')


sf_country <- read_sf('../../metadata/中国_省.geojson')
sf_province <- read_sf('../../metadata/云南省_省.geojson')
sf_city <- read_sf('../../metadata/云南省_市.geojson')

# 定义常用的阿尔伯特投影
albers = "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# 将数据转换为albers投影
sf_country.albers  <- st_transform(sf_country, crs = st_crs(albers))
sf_province.albers <- st_transform(sf_province, crs = st_crs(albers))
sf_city.albers     <- st_transform(sf_city, crs = st_crs(albers))

sampling_point_all <- st_as_sf(cases_data_2019_2022, coords = c("lng", "lat"), crs = st_crs(4326))
sampling_point_all <- st_transform(sampling_point_all, crs = st_crs(albers))

sampling_point_all <- sampling_point_all %>%
  st_join(sf_province.albers, join = st_within, left = FALSE)

sf_china.albers <- sf_country.albers
sf_country.albers <- sf_country.albers %>% dplyr::filter(name != "云南省")

# 中国外边界轮廓
geo_sf.china_outline <- sf_china.albers %>% st_union()



points_df <- as.data.frame(st_coordinates(sampling_point_all))

# 设定带宽参数，采用 5000 米（约对应 10km 的有效范围 2σ）
bw_x <- 50000
bw_y <- 50000

# 扩展边界：上下左右各扩展 2σ（带宽）
expand_x <- range(points_df$X) + c(-2 * bw_x, 2 * bw_x)
expand_y <- range(points_df$Y) + c(-2 * bw_y, 2 * bw_y)

# 使用 kde2d 计算密度，n 控制网格分辨率，使用扩展后的边界
kde_result <- MASS::kde2d(
  x = points_df$X,
  y = points_df$Y,
  n = 1000,
  lims = c(expand_x, expand_y),  # 自定义边界范围
  h = c(bw_x, bw_y)
)

# 将计算结果转换为数据框
density_df <- expand.grid(x = kde_result$x, y = kde_result$y)
density_df$density <- as.vector(kde_result$z)

# 将单位密度转换为 10km×10km 内的样本数
# 注意：面积 = 10,0000 * 10,0000 = 1e10
density_df$sample_count <- density_df$density * 1e10

# 可选：设定阈值，低于一定样本数（如 <1 个样本）区域设为 NA，使之透明
# threshold <- 0.001
# density_df$sample_count[density_df$sample_count < threshold] <- NA

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world.albers <- st_transform(world, crs = st_crs(albers))

# 创建矩形框（原始 CRS = WGS84）
box_proj <- st_as_sfc(st_bbox(c(xmin = 97, xmax = 106, ymin = 21, ymax = 29.1), crs = 4326))
# 投影为地图当前的 CRS（例如 aeqd）
box_proj <- st_transform(box_proj, crs = albers)

ggplot() +
  geom_sf(data = sf_country.albers, fill = 'grey95', linewidth = 0.25, color = "white") +
  geom_sf(data = sf_province.albers, fill = "#ced6e0", color = "black", linewidth = 0.25) + 
  geom_sf(data = geo_sf.china_outline, fill = NA, linewidth = 0.1, color = "black") + 
  geom_sf(data = box_proj, fill = NA, color = "black", linewidth = 0.25, linetype = "dashed") + 
  theme_classic(base_size = 12)

# ggsave('./plot/check_sampling_bias/china_map_background.pdf', width = 3, height = 5)

# 绘制密度图（显示为 10km×10km 内的样本数）
ggplot() + 
  geom_sf(data = sf_province.albers, fill = '#ffffff', color = "black", linewidth = 0.25) + 
  geom_raster(data = density_df, aes(x = x, y = y, fill = sample_count), interpolate = TRUE) + 
  geom_sf(data = world.albers %>% dplyr::filter(sovereignt != "China"), fill = 'white', color = '#ffffff', linewidth = 1) + 
  geom_sf(data = sf_country.albers, fill = "#ffffff", color = "black", linewidth = 0.25) + 
  geom_sf(data = sf_province.albers, fill = NA, color = "black", linewidth = 0.25) + 
  geom_sf(data = sf_city.albers, fill = NA, color = "black", linewidth = 0.1) + 
  # geom_sf(data = sampling_point_all[!is.na(sampling_point_all$tip), ], shape = 21, fill = '#12507b', color = '#ffffff', alpha = 1, size = 2, stroke = 0.5) + 
  scale_fill_gradientn(values = c(0,          0.001,     0.01,      0.1,       0.2,       0.4,       0.8,       1.6), 
                       # colors = c("#ffffff", "#f5f2e9", "#ebe3c7", "#edd3a1", "#f8b862", "#f08300", "#ec6800", "#c9171e"),
                       # colors = c("#ffffff", "#e5f5f5", "#c9e8e8", "#a9dada", "#7fcaca", "#57b5b5", "#389e9e", "#1b8686"),
                       # colors = c("#ffffff", "#e5f1f9", "#c9e2f3", "#a9d0eb", "#7fb9e0", "#57a3d4", "#106898", "#001d52"), 
                       colors = c("grey100", "grey90", "grey80", "grey70", "grey60", "grey50", "grey25", "grey0"), 
                       na.value = NA) + 
  labs(fill = "Total cases\ndensity") + 
  ggnewscale::new_scale_fill() + 
  geom_sf(data = sampling_point_all[!is.na(sampling_point_all$tip), ], 
          aes(fill = group_2, color = group_2), 
          shape = 21, color = '#ffffff', alpha = 1, size = 2.5, stroke = 0.5) + 
  scale_fill_manual(name = "Lineage", values = c("1.1" = '#7EBB8F', "1.2" = '#9DD0CF', "1.3" = '#886A9E')) + 
  labs(fill = "Lineage") + 
  labs(x = "", y = "") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid = element_blank()) + 
  coord_sf(
    crs = albers,   # 投影
    xlim = c(97, 106), # 经纬度
    ylim = c(21, 29.1),
    default_crs = sf::st_crs(4326), # 指示 xlim / ylim 是 EPSG:4326 坐标
    expand = FALSE
  )

# ggsave('./plot/check_sampling_bias/map~sequenced_samples_col2.pdf', width = 5, height = 6)
