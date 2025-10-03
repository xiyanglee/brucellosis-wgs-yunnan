library(tidyverse)
library(geosphere)  # 计算 Haversine 距离
library(Biostrings) # 读取 fasta 文件
library(vegan)      # 计算 Jaccard 距离
library(ape)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v2.RData')



# ------------------------------
# 读取 ERA5 并计算年均值
# ------------------------------

era5_data_monthly <- read.csv('./metadata/sequenced_sample_location_data_with_ecmwf_era5.csv')

colnames(era5_data_monthly)
str(era5_data_monthly)

era5_data_monthly$date  <- as.Date(era5_data_monthly$date)
era5_data_monthly$year  <- year(era5_data_monthly$date)
era5_data_monthly$month <- month(era5_data_monthly$date)
era5_data_monthly$date  <- NULL

era5_data_yearly <- era5_data_monthly %>% 
  group_by(site, year) %>% 
  summarize(
    avg_dewpoint_temperature_2m = mean(dewpoint_temperature_2m), 
    avg_surface_pressure        = mean(surface_pressure), 
    avg_temperature_2m          = mean(temperature_2m), 
    avg_total_precipitation_sum = mean(total_precipitation_sum), 
    .groups = 'drop'
  )

sequenced_cases_data <- cases_data_2019_2022[!is.na(cases_data_2019_2022$tip), ]
sequenced_cases_data$year  <- year(sequenced_cases_data$onset_date)
sequenced_cases_data$month <- month(sequenced_cases_data$onset_date)

sequenced_cases_data <- left_join(sequenced_cases_data, era5_data_monthly, by = c('ID' = 'site', 'year', 'month'))
sequenced_cases_data <- left_join(sequenced_cases_data, era5_data_yearly, by = c('ID' = 'site', 'year'))
rm(era5_data_monthly, era5_data_yearly)

plot(sequenced_cases_data$avg_temperature_2m, sequenced_cases_data$temperature_2m)
plot(sequenced_cases_data$avg_total_precipitation_sum, sequenced_cases_data$total_precipitation_sum)
# surface_pressure 相关性极高，无需区分年度、月度
plot(sequenced_cases_data$avg_surface_pressure, sequenced_cases_data$surface_pressure)
plot(sequenced_cases_data$avg_dewpoint_temperature_2m, sequenced_cases_data$dewpoint_temperature_2m)

# ------------------------------
# 计算测序样本两两之间的各变量差值
# ------------------------------

# 生成所有不重复的 i, j 组合
# combn(seq_len(n), 2) 返回一个 2 x (n*(n-1)/2) 的矩阵
# 每一列对应一对样本索引 (i, j)，确保了不重复且 i < j
# t() 将其转置后得到 n/2 行、2 列，更直观
n         <- nrow(sequenced_cases_data)
pairs_mat <- t(combn(seq_len(n), 2))  # 变为 nC2 行、2 列

# 从 pairs_mat 中分离 i, j
i_idx <- pairs_mat[, 1]
j_idx <- pairs_mat[, 2]

## 向量化距离计算
## distHaversine() 支持将坐标输入为两列矩阵（每行对应一组坐标），会返回一个向量，向量长度与行数（即对数）匹配

# 2. 逐列向量化计算距离与差值
# 2.1 计算地理距离（米）再转换为 km
dist_m <- distHaversine(
  cbind(sequenced_cases_data$lng[i_idx], sequenced_cases_data$lat[i_idx]),
  cbind(sequenced_cases_data$lng[j_idx], sequenced_cases_data$lat[j_idx])
)
dist_km <- dist_m / 1000

# 2.2 计算时间间隔（天），取绝对值
time_interval <- abs(
  as.numeric(
    difftime(
      sequenced_cases_data$onset_date[j_idx],
      sequenced_cases_data$onset_date[i_idx],
      units = "days"
    )
  )
)

# 2.3 计算各项数值的绝对差
monthly_dewpoint_diff <- abs(
  sequenced_cases_data$dewpoint_temperature_2m[i_idx] -
    sequenced_cases_data$dewpoint_temperature_2m[j_idx]
)

monthly_surface_pressure_diff <- abs(
  sequenced_cases_data$surface_pressure[i_idx] -
    sequenced_cases_data$surface_pressure[j_idx]
)

monthly_temp_2m_diff <- abs(
  sequenced_cases_data$temperature_2m[i_idx] -
    sequenced_cases_data$temperature_2m[j_idx]
)

monthly_precip_diff <- abs(
  sequenced_cases_data$total_precipitation_sum[i_idx] -
    sequenced_cases_data$total_precipitation_sum[j_idx]
)

yearly_dewpoint_diff <- abs(
  sequenced_cases_data$avg_dewpoint_temperature_2m[i_idx] -
    sequenced_cases_data$avg_dewpoint_temperature_2m[j_idx]
)

yearly_surface_pressure_diff <- abs(
  sequenced_cases_data$avg_surface_pressure[i_idx] -
    sequenced_cases_data$avg_surface_pressure[j_idx]
)

yearly_temp_diff <- abs(
  sequenced_cases_data$avg_temperature_2m[i_idx] -
    sequenced_cases_data$avg_temperature_2m[j_idx]
)

yearly_precip_diff <- abs(
  sequenced_cases_data$avg_total_precipitation_sum[i_idx] -
    sequenced_cases_data$avg_total_precipitation_sum[j_idx]
)

# 3. 将全部计算结果一次性组合成新的数据框
pairwise_sequenced_cases_data <- data.frame(
  tip_pair  = paste0(sequenced_cases_data$tip[i_idx], "_", sequenced_cases_data$tip[j_idx]),
  tip_1     = sequenced_cases_data$tip[i_idx],
  tip_2     = sequenced_cases_data$tip[j_idx],
  lineage_1 = sequenced_cases_data$group_2[i_idx],
  lineage_2 = sequenced_cases_data$group_2[j_idx],
  
  distance_km  = dist_km,
  time_interval = time_interval,
  
  monthly_dewpoint_temperature_2m_difference   = monthly_dewpoint_diff,
  monthly_surface_pressure_2m_difference       = monthly_surface_pressure_diff,
  monthly_temperature_2m_2m_difference         = monthly_temp_2m_diff,
  monthly_total_precipitation_sum_difference   = monthly_precip_diff,
  
  yearly_dewpoint_temperature_2m_difference    = yearly_dewpoint_diff,
  yearly_surface_pressure_2m_difference        = yearly_surface_pressure_diff, 
  yearly_temperature_2m_difference             = yearly_temp_diff,
  yearly_total_precipitation_sum_difference    = yearly_precip_diff,
  
  stringsAsFactors = FALSE
)

rm(n, pairs_mat, i_idx, j_idx,
   dist_m, dist_km, time_interval,
   monthly_dewpoint_diff, monthly_surface_pressure_diff,
   monthly_temp_2m_diff, monthly_precip_diff,
   yearly_dewpoint_diff, yearly_temp_diff, yearly_precip_diff)

# ------------------------------
# 计算两两样本 平均核苷酸相似性（ANI，Average Nucleotide Identity）
# ------------------------------

# 1. 读取序列并提取样本名与序列
fasta_file   <- readDNAStringSet("./data/snippy_contig_gubbins/clean.core.aln")
sample_names <- names(fasta_file)                 # 与 pairwise_sequenced_cases_data 中的 tip 对应
sequences    <- as.character(fasta_file)          # 将 DNAString 转换为普通字符向量，每个元素是一整条序列

# 2. 建立一个命名向量，方便通过样本名快速找到序列
seq_dict <- setNames(sequences, sample_names)

# 3. 逐行计算 ANI
#   - 这里直接用 apply(...) 对数据框按行处理
#   - 对于每一行，通过 row["tip_1"]、row["tip_2"] 找到各自序列
#   - strsplit(...) 切成字符向量后用 == 比较
#   - mean(...) 即为碱基相同的比例，视为 ANI
pairwise_sequenced_cases_data$ANI <- apply(pairwise_sequenced_cases_data, 1, function(row) {
  seq1 <- seq_dict[[ row["tip_1"] ]]
  seq2 <- seq_dict[[ row["tip_2"] ]]
  
  # 将字符串拆为碱基向量
  seq1_vec <- strsplit(seq1, "")[[1]]
  seq2_vec <- strsplit(seq2, "")[[1]]
  
  # 计算平均相似度
  mean(seq1_vec == seq2_vec)
})

rm(fasta_file, sample_names, sequences, seq_dict)

# ------------------------------
# 基因含量相似度（Gene content similarity），即 Jaccard 距离
# ------------------------------

# datmat 是行名为样本 tip、列名为基因名、值为 0/1 的矩阵
datmat <- read.table('./data/roary_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat <- t(datmat) %>% as.data.frame()

colnames(datmat) <- datmat[1, ]
datmat <- datmat[-1, ]

selected_indices <- grepl("^20", datmat$Gene)  # 筛选样本名以 "20" 开头的行
datmat <- datmat[selected_indices, ]

rownames_datamat <- datmat$Gene
datmat$Gene <- NULL
datmat <- datmat %>% reframe(across(everything(), as.integer))
rownames(datmat) <- rownames_datamat

# 2. 计算 Jaccard 距离矩阵，然后转成相似度矩阵
distmat <- vegdist(datmat, method = "jaccard", binary = TRUE)
distmat <- as.matrix(distmat)
gene_sim <- 1 - distmat  # 将距离转成相似度

# 3. 按行将每对 (tip_1, tip_2) 的相似度写回 pairwise_sequenced_cases_data
#   这里用 mapply 来避免显式的循环。
pairwise_sequenced_cases_data$gene_content_similarity <- mapply(
  FUN = function(a, b) gene_sim[a, b],
  pairwise_sequenced_cases_data$tip_1,
  pairwise_sequenced_cases_data$tip_2
)

rm(distmat, gene_sim, rownames_datamat, selected_indices)

# ------------------------------
# Phylogenetic similarity 计算
# ------------------------------

# (1) 计算所有 tips 的距离矩阵
dist_matrix <- cophenetic.phylo(raxml_tree)

# (2) 转化为相似度矩阵
max_dist           <- max(dist_matrix)
phylo_sim_matrix   <- 1 - dist_matrix / max_dist

# (3) 利用 mapply 索引并填充到新列
pairwise_sequenced_cases_data$phylo_similarity <- mapply(
  function(a, b) phylo_sim_matrix[a, b],
  pairwise_sequenced_cases_data$tip_1,
  pairwise_sequenced_cases_data$tip_2
)

# 如需删除中间变量
rm(dist_matrix, max_dist, phylo_sim_matrix)



# save.image('./r_image/compare_lineage_yunnan_v2.RData')


