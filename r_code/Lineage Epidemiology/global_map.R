# install.packages('scatterpie', lib = "/public/r_share_library/4.1")
library(dplyr)
library(ggplot2)
library(scatterpie)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

setwd("/Users/xiyangli/Lab/Project/Brucella_WGS_Yunnan_2024/pipelines/denovo_ssembly")

load('./r_image/compare_lineage_yunnan_v2.RData')



# > table(metadata_tree %>% select(area, group_2))
# group_2
# area              1.1 1.2 1.3 2.1 2.2 2.3 3.1 3.2
# Central Asia      2   0   0   0   0   0   0   0
# Central Europe    3   9   0   0   0   0   0   5
# East Asia       140   3   1   0   0   1   0   0
# Eastern Africa    1   1   0   5   0   0   0   0
# Eastern Europe    3   0   0   0   0   0   0   0
# North Africa      0   0   0   0   0   0   0   4
# North America     0   1   0   0   0   1   1   1
# North Asia       23   0   0   0   0   0   0   0
# North Europe      5  13   0   6   0   1   0   0
# South America     0   0   0   0   0   3   0   0
# South Asia       42   0   0   0   0   1   0   0
# Southeast Asia   14   0   0   0   0   0   0   0
# Southern Africa   0   0   0   3  12   1   0   0
# Southern Europe   0   8   0   0   0   2  90  22
# Western Africa    0   0   0   1   0   0   0   0
# Western Asia      8 256   1   1   0   0   0   0
# Western Europe    2   1   0   2   0   1   7   0
# Yunnan           79   5  19   0   0   0   0   0

# 获取世界地图
world <- ne_countries(scale = "medium", returnclass = "sf")
world_coastline <- ne_coastline(scale = "medium", returnclass = "sf")

# 读取中国地图
china_sf <- read_sf('./source/tianditu/中国_省.geojson') %>% dplyr::filter(name != '境界线')

# 自定义每个区域的大致经纬度中心点
area_coords <- tibble::tribble(
  ~area,               ~lon,   ~lat,
  "Yunnan",            110,    25,
  "East Asia",         130,    35,
  "Southeast Asia",    105,    10,
  "South Asia",         80,    22,
  "Central Asia",       70,    45,
  "Western Asia",       45,    33,
  "North Asia",        100,    60,
  "Eastern Europe",     30,    55,
  "Central Europe",     15,    50,
  "Western Europe",      5,    48,
  "Southern Europe",    15,    41,
  "North Europe",       20,    60,
  "North America",    -100,    45,
  "South America",     -60,   -15,
  "North Africa",       10,    30,
  "Western Africa",     -5,    10,
  "Eastern Africa",     40,     0,
  "Southern Africa",    25,   -30
)

# 汇总饼图数据
pie_data <- metadata_tree %>%
  group_by(area, group_2) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = group_2, values_from = n, values_fill = 0)

# 计算每个 area 的总样本数，用作饼图直径比例
pie_data <- pie_data %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  ungroup()

# 合并经纬度信息
pie_data <- left_join(pie_data, area_coords, by = "area")

# 统计每个国家的样本数，去除 NA
country_counts <- metadata_tree %>% 
  mutate(geographic_country = case_match(
    geographic_country, 
    "Yunnan" ~ "China", 
    .default = geographic_country
  )) %>% 
  filter(!is.na(geographic_country)) %>%
  count(geographic_country, name = "sample_n")

# 合并到世界地图数据中，匹配列名可根据 world 的实际国家名列（如 'name'）调整
world <- world %>% 
  left_join(country_counts, by = c("name" = "geographic_country")) 

# 添加分组列
world <- world %>%
  mutate(sample_cat = case_when(
    is.na(sample_n)      ~ "0",
    sample_n <= 10       ~ "0–10",
    sample_n <= 50       ~ "10–50",
    sample_n <= 100      ~ "50–100",
    sample_n > 100       ~ ">100"
  ))


# 创建边界框并投影为 Equal Earth
grat1 <- st_graticule(lat = seq(-90, 90, by = 179),
                     lon = seq(-180, 180, by = 360),
                     crs = "+proj=eqearth")
grat2 <- st_graticule(lat = seq(90, -90, by = -179),
                     lon = seq(-180, 180, by = 360),
                     crs = "+proj=eqearth")

ggplot() + 
  geom_sf(data = world_coastline, fill = "gray95", color = 'black', linewidth = 0.25) + 
  geom_sf(data = world, aes(fill = sample_cat), color = NA, size = 0.2, show.legend = TRUE) + 
  geom_sf(data = china_sf, fill = "#444444", color = NA) + 
  scale_fill_manual(
    values = c(
      "0–10" = "#dddddd",       # 浅灰
      "10–50" = "#bbbbbb",      # 中灰
      "50–100" = "#888888",     # 深灰
      ">100" = "#444444"        # 很深灰
    ),
    na.value = "white",
    name = "Sample count",
    drop = FALSE
  ) + 
  geom_sf(data = grat1, color = "black", size = 0.3, fill = NA) +  # 投影边缘来自经纬网
  geom_sf(data = grat2, color = "black", size = 0.3, fill = NA) +  # 投影边缘来自经纬网
  theme_void() + 
  coord_sf(crs = "+proj=eqearth")

# ggsave('./plot/global_map/global_map_background_eqearth.pdf', width = 10)

ggplot() + 
  geom_scatterpie(
    aes(x = lon, y = lat, group = area, r = log1p(sqrt(total)) * 4),
    data = pie_data,
    cols = setdiff(names(pie_data), c("area", "lon", "lat", "total")),  # 分类列
    color = "black",  linewidth = 0.1,   # 饼图边界
    alpha = 1
  ) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
                               "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
                               "3.1" = "#c9a7b7", "3.2" = "#EBCACB"), 
                    na.value = "transparent") + 
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) + 
  geom_scatterpie_legend(radius = log1p(sqrt(100)) * 4, labeller = function(.) 100, x = -120, y = 0) + 
  geom_scatterpie_legend(radius = log1p(sqrt(50)) * 4, labeller = function(.) 50, x = -120, y = -20) + 
  geom_scatterpie_legend(radius = log1p(sqrt(10)) * 4, labeller = function(.) 10, x = -120, y = -40) + 
  coord_fixed()

# ggsave('./plot/global_map/global_map_pie_only.pdf', width = 10)



ggplot() + 
  geom_sf(data = world, aes(fill = sample_cat), color = NA, size = 0.1, show.legend = FALSE) + 
  geom_sf(data = world_coastline, fill = "gray95", color = 'black', linewidth = 0.15) +
  geom_sf(data = china_sf, fill = "#444444", color = NA) + 
  scale_fill_manual(
    values = c(
      "0–10" = "#dddddd",
      "10–50" = "#bbbbbb",
      "50–100" = "#888888",
      ">100" = "#444444"
    ),
    na.value = "#ffffff",
    name = "Sample count",
    drop = FALSE
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_scatterpie(
    aes(x = lon, y = lat, group = area, r = log1p(sqrt(total)) * 4),
    data = pie_data,
    cols = setdiff(names(pie_data), c("area", "lon", "lat", "total")),
    color = "black", linewidth = 0.25, alpha = 1
  ) + 
  scale_fill_manual(
    values = c(
      "1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
      "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
      "3.1" = "#c9a7b7", "3.2" = "#EBCACB"
    ), 
    na.value = "transparent"
  ) +
  coord_sf(
    # crs = paste0("+proj=ortho +lon_0=60 +lat_0=15 +datum=WGS84"),
    crs = paste0("+proj=ortho +lon_0=-100 +lat_0=15 +datum=WGS84"),
    clip = "on"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA)
  )

# ggsave('./plot/global_map/global_map_orthographic_projection_lon_60.pdf', width = 8, height = 8)
# ggsave('./plot/global_map/global_map_orthographic_projection_lon_-100.pdf', width = 8, height = 8)




# 定义信息熵函数
calc_entropy <- function(x) {
  freq <- table(x) / length(x)
  -sum(freq * log2(freq))
}

entropy_df <- metadata_tree %>%
  group_by(area) %>%
  summarise(
    n_sample = n(),
    entropy = calc_entropy(group_2)
  ) %>%
  ungroup()

lm_fit <- lm(entropy ~ n_sample, data = entropy_df)
newdata <- data.frame(n_sample = seq(min(entropy_df$n_sample) - 0.2, max(entropy_df$n_sample) + 50, length.out = 100))
pred <- predict(lm_fit, newdata = newdata, se.fit = TRUE)
newdata$fit <- pred$fit
newdata$lwr <- pred$fit - 1.96 * pred$se.fit
newdata$upr <- pred$fit + 1.96 * pred$se.fit

summary(lm_fit)

# Call:
#   lm(formula = entropy ~ n_sample, data = entropy_df)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -0.7795 -0.7723 -0.2741  0.5679  1.4201 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  0.7804082  0.2239667   3.484  0.00284 **
# n_sample    -0.0009042  0.0027891  -0.324  0.74976   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8124 on 17 degrees of freedom
# Multiple R-squared:  0.006144,	Adjusted R-squared:  -0.05232 
# F-statistic: 0.1051 on 1 and 17 DF,  p-value: 0.7498

ggplot() + 
  geom_ribbon(data = newdata, aes(x = n_sample, ymin = lwr, ymax = upr), alpha = 0.1, color = "black", linetype = "dashed") + 
  geom_line(data = newdata, aes(x = n_sample, y = fit), linewidth = 1) + 
  geom_point(data = entropy_df, aes(x = n_sample, y = entropy), size = 4, color = "black") + 
  geom_point(data = entropy_df, aes(x = n_sample, y = entropy, color = area), size = 3) + 
  # geom_smooth(method = "lm") + 
  labs(x = "Number of samples", y = "Shannon entropy of group_2") + 
  scale_x_log10() + 
  coord_cartesian(ylim = c(-0.1, 2.1), expand = 0) + 
  theme_classic(base_size = 14) + 
  scale_color_manual(values = c("East Asia" = "#4f6f46", 
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
  annotation_logticks(base = 10, sides = "b") + 
  theme(legend.position = "none")

# ggsave('./plot/global_map/entropy_per_area.pdf', width = 4.5, height = 4.5)

ggplot() + 
  # geom_ribbon(data = newdata, aes(x = n_sample, ymin = lwr, ymax = upr), alpha = 0.1, color = "black", linetype = "dashed") + 
  # geom_line(data = newdata, aes(x = n_sample, y = fit), linewidth = 1) + 
  geom_point(data = entropy_df, aes(x = n_sample, y = entropy), size = 4, color = "black") + 
  geom_point(data = entropy_df, aes(x = n_sample, y = entropy, color = area), size = 3) + 
  # geom_smooth(method = "lm") + 
  labs(x = "Number of samples", y = "Shannon entropy of group_2") + 
  scale_x_log10() + 
  # coord_cartesian(ylim = c(-0.1, 2.1), expand = 0) + 
  theme_classic(base_size = 14) + 
  scale_color_manual(values = c("East Asia" = "#4f6f46", 
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
  annotation_logticks(base = 14, sides = "b") + 
  theme(legend.position = "none")

# ggsave('./plot/global_map/entropy_per_area_type2.pdf', width = 3.5, height = 3.5)

entropy_df %>% 
  na.omit() %>% 
  write.csv("./r_plot_data/figure_1c~entropy_df.csv", row.names = FALSE)
