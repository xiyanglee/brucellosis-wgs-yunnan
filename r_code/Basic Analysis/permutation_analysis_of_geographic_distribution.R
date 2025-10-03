# install.packages('geosphere', lib = '/public/r_share_library/4.1')

library(sf)         # 处理矢量数据
library(ggplot2)    # 创建地图
library(geosphere)  # 计算 Haversine 距离
library(parallel)   # 多线程
library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v2.RData')


sf_province <- read_sf('../../metadata/云南省_省.geojson')
# 定义阿尔伯特投影
albers = "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# 将数据转换为 albers 投影
sf_province.albers  <- st_transform(sf_province, crs = st_crs(albers))

sampling_point_all       <- cases_data_2019_2022
sampling_point_all$lng_1 <- sampling_point_all$lng
sampling_point_all$lat_1 <- sampling_point_all$lat
sampling_point_all       <- st_as_sf(sampling_point_all, coords = c("lng_1", "lat_1"), crs = st_crs(4326))
sampling_point_all       <- st_transform(sampling_point_all, crs = st_crs(albers))
sampling_point_all       <- sampling_point_all %>%
  st_join(sf_province.albers, join = st_within, left = FALSE)
sampling_point_all       <- as.data.frame(sampling_point_all)

sequencing_point         <- sampling_point_all[!is.na(sampling_point_all$tip), ]
sequencing_point         <- as.data.frame(sequencing_point)

# ------------------------------
# 1. 定义计算伪F值的函数
# ------------------------------
compute_pseudoF <- function(df, group){
  # df: 包含 lat 和 lng 列的 data.frame
  # group: 与 df 行数相同的分组因子
  #
  # 计算两两样本之间的 Haversine 距离（单位：公里），并取平方后计算总平方和
  # 再按照组内、组间进行分割，计算伪F值
  
  # 注意：geosphere::distm 要求坐标矩阵顺序为 (lng, lat)
  coords <- as.matrix(df[, c("lng", "lat")])
  dist_mat <- distm(coords, fun = distHaversine)
  dist_sq <- dist_mat^2  # 距离的平方
  
  N <- nrow(df)
  # 总平方和：仅计算上三角（不重复计算）
  SS_total <- sum(dist_sq[upper.tri(dist_sq)])
  
  # 计算各组内部的平方和
  groups <- unique(group)
  SS_within <- 0
  for(g in groups){
    idx <- which(group == g)
    if(length(idx) > 1){
      sub_mat <- dist_sq[idx, idx]
      SS_within <- SS_within + sum(sub_mat[upper.tri(sub_mat)])
    }
  }
  # 组间平方和
  SS_between <- SS_total - SS_within
  
  a <- length(groups)  # 组数（应为 2）
  pseudoF <- (SS_between/(a - 1)) / (SS_within/(N - a))
  return(pseudoF)
}

# ------------------------------
# 2. 数据准备
# ------------------------------
# sampling_point_all 为所有病例的 data.frame，包含 lat, lng, tip 列
# tip 不为 NA 的为已测序样本，tip 为 NA 的为未测序样本
# 构建观测分组
group_obs <- ifelse(is.na(sampling_point_all$tip), "nonsequenced", "sequenced")
group_obs <- factor(group_obs, levels = c("sequenced", "nonsequenced"))

# ------------------------------
# 3. 计算观测 pseudoF 值（原始分组）
# ------------------------------
observed_pseudoF <- compute_pseudoF(sampling_point_all, group_obs)
cat("Observed pseudoF:", observed_pseudoF, "\n")

# 记录原始已测序样本数量
n_seq <- sum(!is.na(sampling_point_all$tip))

# ------------------------------
# 4. 置换检验（多线程实现）
# ------------------------------
n_perm <- 999
set.seed(123)  # 固定随机种子保证结果可重复

# 使用 mclapply 进行多线程计算（注意：Windows 用户可以考虑使用 parLapply 创建集群）
perm_results <- mclapply(1:n_perm, function(i) {
  # 随机抽取 n_seq 个病例作为置换“sequenced”组，其余划为对照组
  perm_indices <- sample(1:nrow(sampling_point_all), n_seq, replace = FALSE)
  group_perm <- rep("control", nrow(sampling_point_all))
  group_perm[perm_indices] <- "sequenced"
  group_perm <- factor(group_perm, levels = c("sequenced", "control"))
  
  compute_pseudoF(sampling_point_all, group_perm)
}, mc.cores = detectCores())

perm_pseudoF <- unlist(perm_results)

# 计算 p-value
p_value <- (sum(perm_pseudoF >= observed_pseudoF) + 1) / (n_perm + 1)
cat("Permutation test p-value:", p_value, "\n")

# 结果解释：
# 如果 p-value > 0.05，则说明在至少 5% 的随机置换中计算得到的伪F值不低于观测值，
# 表明两组（测序与非测序）的地理分布无显著差异，即测序样本在地理上具有代表性。


ggplot(data.frame(PseudoF = perm_pseudoF), aes(x = PseudoF)) +
  geom_histogram(bins = 50, color = "#6e9bc5", fill = "#aed0ee") +
  geom_vline(xintercept = observed_pseudoF, color = "#c12c1f", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Permutation Pseudo-F Values",
       x = "Pseudo-F Value", y = "Frequency") +
  theme_classic(base_size = 14)

# ggsave('./plot/check_sampling_bias/hist~perm_pseudo_F.pdf', width = 5, height = 3)

ggplot(data.frame(PseudoF = perm_pseudoF), aes(x = PseudoF)) +
  stat_ecdf(geom = "step", color = "#6e9bc5", linewidth = 1) +
  geom_vline(xintercept = observed_pseudoF, color = "#c12c1f", linetype = "dashed", linewidth = 1) +
  labs(title = "ECDF of Permutation Pseudo-F Values",
       x = "Pseudo-F Value", y = "Cumulative Proportion") +
  theme_bw(base_size = 14)

# ggsave('./plot/check_sampling_bias/line~ECDF_of_pseudo_F.pdf', width = 4, height = 3.5)

data.frame(
  Mean_PseudoF = mean(perm_pseudoF),
  Median_PseudoF = median(perm_pseudoF),
  SD_PseudoF = sd(perm_pseudoF),
  Observed_PseudoF = observed_pseudoF,
  p_value = p_value
)

# Mean_PseudoF Median_PseudoF SD_PseudoF Observed_PseudoF p_value
#     217.3081       216.7641   5.629489         213.8842   0.713

