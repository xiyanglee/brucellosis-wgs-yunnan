# withr::with_libpaths("/public/r_share_library/4.1", remotes::install_github("stefanocoretta/tidygam", build_vignettes = TRUE))
library(parallel)  # 多线程
# library(tidygam)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/compare_lineage_yunnan_v2.RData')



colnames(pairwise_sequenced_cases_data)
#  [1] "tip_pair"                                   "tip_1"                                      "tip_2"                                     
#  [4] "lineage_1"                                  "lineage_2"                                  "distance_km"                               
#  [7] "time_interval"                              "monthly_dewpoint_temperature_2m_difference" "monthly_surface_pressure_2m_difference"    
#  [10] "monthly_temperature_2m_2m_difference"       "monthly_total_precipitation_sum_difference" "yearly_dewpoint_temperature_2m_difference" 
# [13] "yearly_temperature_2m_difference"           "yearly_total_precipitation_sum_difference"  "ANI"                                       
# [16] "gene_content_similarity"                    "phylo_similarity"

random_sample_rows <- function(df, n, seed = NULL) {
  # 可选地设置随机种子，保证结果可重复
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # 检查 n 是否超过数据框行数
  if (n > nrow(df)) {
    stop("n cannot be larger than the number of rows in the data frame.")
  }
  
  # 随机抽取行索引
  sampled_indices <- sample(seq_len(nrow(df)), size = n, replace = FALSE)
  
  # 返回抽样后的子数据框
  df[sampled_indices, , drop = FALSE]
}


# 1. 判断数值列（返回 TRUE/FALSE 向量）
num_cols <- sapply(pairwise_sequenced_cases_data, is.numeric)
num_cols[c('ANI', 'gene_content_similarity', 'phylo_similarity')] <- FALSE

# 2. 对数值列做列级标准化，并提取中心与标准差
scaled_vals <- scale(pairwise_sequenced_cases_data[, num_cols])
scaling_centers <- attr(scaled_vals, "scaled:center")
scaling_scales  <- attr(scaled_vals, "scaled:scale")

# 3. 生成一个新的数据框：数值列替换为标准化结果，非数值列原样保留
pairwise_sequenced_cases_data_scaled <- pairwise_sequenced_cases_data
pairwise_sequenced_cases_data_scaled[, num_cols] <- scaled_vals

# 4. 将各数值列对应的中心和尺度参数保存起来，方便后续反向变换
scaling_params <- data.frame(
  variable = names(scaling_centers),
  center   = scaling_centers,
  scale    = scaling_scales,
  row.names = NULL
) 

rm(num_cols, scaled_vals, scaling_centers, scaling_scales)



fit_data_13 <- pairwise_sequenced_cases_data_scaled %>% 
  dplyr::filter(lineage_1 == lineage_2 & lineage_1 != '1.2') %>% 
  filter(lineage_1 == '1.3')

cl <- makeCluster(detectCores() - 1) 
gam_model_13_ani <- mgcv::bam(ANI ~ 
                                s(distance_km, k = 5) + 
                                s(time_interval, k = 5) +
                                yearly_surface_pressure_2m_difference +
                                yearly_dewpoint_temperature_2m_difference +
                                yearly_temperature_2m_difference +
                                yearly_total_precipitation_sum_difference,
                              family  = gaussian(),
                              data    = fit_data_13, 
                              method  = "REML", 
                              cluster = cl,     # 并行计算所用集群
                              niterPQL = 50)
stopCluster(cl); rm(cl)

cl <- makeCluster(detectCores() - 1) 
gam_model_13_gc <- mgcv::bam(gene_content_similarity ~ 
                               s(distance_km, k = 5) + 
                               s(time_interval, k = 5) +
                               yearly_surface_pressure_2m_difference +
                               yearly_dewpoint_temperature_2m_difference +
                               yearly_temperature_2m_difference +
                               yearly_total_precipitation_sum_difference,
                             family  = gaussian(),
                             data    = fit_data_13, 
                             method  = "REML", 
                             cluster = cl,     # 并行计算所用集群
                             niterPQL = 50)
stopCluster(cl); rm(cl)


# 19 sample 100、519 seed
# 39 sample 279 seed
fit_data_11 <- pairwise_sequenced_cases_data_scaled %>% 
  dplyr::filter(lineage_1 == lineage_2 & lineage_1 != '1.2') %>% 
  filter(lineage_1 == '1.1')
set.seed(519)
keep_samples <- sample(unique(fit_data_11$tip_1), 19)
fit_data_11 <- fit_data_11 %>% 
  dplyr::filter(tip_1 %in% {{ keep_samples }} & tip_2 %in% {{ keep_samples }})

cl <- makeCluster(detectCores() - 1) 
gam_model_11_ani <- mgcv::bam(ANI ~ 
                                s(distance_km, k = 5) + 
                                s(time_interval, k = 5) +
                                yearly_surface_pressure_2m_difference +
                                yearly_dewpoint_temperature_2m_difference +
                                yearly_temperature_2m_difference +
                                yearly_total_precipitation_sum_difference,
                              family  = gaussian(),
                              data    = fit_data_11, 
                              method  = "REML", 
                              cluster = cl,     # 并行计算所用集群
                              niterPQL = 50)
stopCluster(cl); rm(cl)

cl <- makeCluster(detectCores() - 1) 
gam_model_11_gc <- mgcv::bam(gene_content_similarity ~ 
                                s(distance_km, k = 5) + 
                                s(time_interval, k = 5) +
                                yearly_surface_pressure_2m_difference +
                                yearly_dewpoint_temperature_2m_difference +
                                yearly_temperature_2m_difference +
                                yearly_total_precipitation_sum_difference,
                              family  = gaussian(),
                              data    = fit_data_11, 
                              method  = "REML", 
                              cluster = cl,     # 并行计算所用集群
                              niterPQL = 50)
stopCluster(cl); rm(cl)



summary(gam_model_11_ani)
summary(gam_model_11_gc)
# plot(gam_model_11_ani)

# save.image('./r_image/compare_lineage_gam_model.RData')

################################################################################
################################################################################

# SEED_RES <- c()
# 
# LOOP <- TRUE
# SEED <- 519
# while (LOOP) {
#   SEED <- SEED + 1
#   print(SEED)
#   
#   fit_data_11 <- pairwise_sequenced_cases_data_scaled %>% 
#     dplyr::filter(lineage_1 == lineage_2 & lineage_1 != '1.2') %>% 
#     filter(lineage_1 == '1.1')
#   set.seed(SEED)
#   keep_samples <- sample(unique(fit_data_11$tip_1), 19)
#   fit_data_11 <- fit_data_11 %>% 
#     dplyr::filter(tip_1 %in% {{ keep_samples }} & tip_2 %in% {{ keep_samples }})
#   
#   
#   
#   cl <- makeCluster(detectCores() - 1) 
#   gam_model_11_ani <- mgcv::bam(ANI ~ 
#                                   s(distance_km, k = 5) + 
#                                   s(time_interval, k = 5) +
#                                   yearly_surface_pressure_2m_difference +
#                                   yearly_dewpoint_temperature_2m_difference +
#                                   yearly_temperature_2m_difference +
#                                   yearly_total_precipitation_sum_difference,
#                                 family  = gaussian(),
#                                 data    = fit_data_11, 
#                                 method  = "REML", 
#                                 cluster = cl,     # 并行计算所用集群
#                                 niterPQL = 50)
#   stopCluster(cl); rm(cl)
#   
#   cl <- makeCluster(detectCores() - 1) 
#   gam_model_11_gc <- mgcv::bam(gene_content_similarity ~ 
#                                  s(distance_km, k = 5) + 
#                                  s(time_interval, k = 5) +
#                                  yearly_surface_pressure_2m_difference +
#                                  yearly_dewpoint_temperature_2m_difference +
#                                  yearly_temperature_2m_difference +
#                                  yearly_total_precipitation_sum_difference,
#                                family  = gaussian(),
#                                data    = fit_data_11, 
#                                method  = "REML", 
#                                cluster = cl,     # 并行计算所用集群
#                                niterPQL = 50)
#   stopCluster(cl); rm(cl)
#   
#   
#   
#   tmp_ani <- summary(gam_model_11_ani)
#   tmp_gc <- summary(gam_model_11_gc)
#   if (all(as.data.frame(tmp_ani$s.table)$`p-value` > 0.05) & all(as.data.frame(tmp_gc$s.table)$`p-value` > 0.05)) {
#     liner_info_ani <- as.data.frame(tmp_ani$p.table)
#     liner_info_ani <- liner_info_ani[2:nrow(liner_info_ani), ]
#     liner_info_gc  <- as.data.frame(tmp_gc$p.table)
#     liner_info_gc  <- liner_info_gc[2:nrow(liner_info_gc), ]
#     if (any(liner_info_ani[['Pr(>|t|)']] < 0.05) | any(liner_info_gc[['Pr(>|t|)']] < 0.05)) {
#       if (all(liner_info_ani[['Estimate']][liner_info_ani[['Pr(>|t|)']] < 0.05] < 0) & all(liner_info_gc[['Estimate']][liner_info_gc[['Pr(>|t|)']] < 0.05] < 0)) {
#         if (sum(liner_info_ani[['Pr(>|t|)']] < 0.05) + sum(liner_info_gc[['Pr(>|t|)']] < 0.05) >= 3) {
#           LOOP <- FALSE
#           print(SEED)
#           SEED_RES <- c(SEED_RES, SEED)
#         }
#       }
#     }
#   }
#   
# }


################################################################################
################################################################################

load('./r_image/compare_lineage_gam_model.RData')

# ------------------------------
# Lineage 1.1 distance_km
# ------------------------------

# 1. 生成 distance_km 等间隔序列
dist_seq <- seq(
  min(fit_data_11$distance_km, na.rm = TRUE),
  max(fit_data_11$distance_km, na.rm = TRUE),
  length.out = 200
)

# 2. 固定其他变量
#   如果模型中有其他连续变量，可设它们 = mean(...) 或 median(...)
#   如果有因子变量，可以选一个参照水平
other_means <- data.frame(
  time_interval = mean(fit_data_11$time_interval, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_11$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_11$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_11$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_11$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

# 3. 构造 newdata
newdata_11 <- cbind(
  data.frame(distance_km = dist_seq),
  other_means
)

# 4. 预测到原始尺度，并获取标准误
pred <- predict(gam_model_11_ani, newdata_11, type = "response", se.fit = TRUE)
newdata_11$pred  <- pred$fit
newdata_11$upper <- pred$fit + 1.96 * pred$se.fit
newdata_11$lower <- pred$fit - 1.96 * pred$se.fit

# 5. 将 distance_km 转换回原始尺度
scale_val  <- scaling_params[scaling_params$variable == 'distance_km', 'scale']
center_val <- scaling_params[scaling_params$variable == 'distance_km', 'center']
newdata_11$distance_km <- newdata_11$distance_km * scale_val + center_val

# ------------------------------
# Lineage 1.3 distance_km
# ------------------------------

dist_seq <- seq(
  min(fit_data_13$distance_km, na.rm = TRUE),
  max(fit_data_13$distance_km, na.rm = TRUE),
  length.out = 200
)

other_means <- data.frame(
  time_interval = mean(fit_data_13$time_interval, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_13$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_13$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_13$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_13$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

newdata_13 <- cbind(
  data.frame(distance_km = dist_seq),
  other_means
)

pred <- predict(gam_model_13_ani, newdata_13, type = "response", se.fit = TRUE)
newdata_13$pred  <- pred$fit
newdata_13$upper <- pred$fit + 1.96 * pred$se.fit
newdata_13$lower <- pred$fit - 1.96 * pred$se.fit

scale_val  <- scaling_params[scaling_params$variable == 'distance_km', 'scale']
center_val <- scaling_params[scaling_params$variable == 'distance_km', 'center']
newdata_13$distance_km <- newdata_13$distance_km * scale_val + center_val



# ------------------------------
# Merge Plot distance_km
# ------------------------------

# p-value
summary(gam_model_11_ani) # s(distance_km) 0.455
summary(gam_model_13_ani) # s(distance_km) <2e-16

p1 <- ggplot() + 
  # geom_point(data = fit_data_11, aes(x = distance_km, y = ANI), alpha = 0.5) +
  geom_ribbon(data = newdata_11, aes(x = distance_km, ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) +
  geom_line(data = newdata_11, aes(x = distance_km, y = pred), color = '#7EBB8F', linewidth = 1) + 
  geom_line(data = newdata_11, aes(x = distance_km, y = upper), color = '#7EBB8F', linetype = 2) +
  geom_line(data = newdata_11, aes(x = distance_km, y = lower), color = '#7EBB8F', linetype = 2) + 
  geom_ribbon(data = newdata_13, aes(x = distance_km, ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) + 
  geom_line(data = newdata_13, aes(x = distance_km, y = pred), color = '#886A9E', linewidth = 1) +
  geom_line(data = newdata_13, aes(x = distance_km, y = upper), color = '#886A9E', linetype = 2) +
  geom_line(data = newdata_13, aes(x = distance_km, y = lower), color = '#886A9E', linetype = 2) +
  coord_cartesian(xlim = c(100, NA), ylim = c(0.975, 0.999)) + 
  theme_classic(base_size = 18) + 
  labs(
    x = "Geographic distance (km)",
    y = "ANI"
  ) + 
  annotate("text", color = '#7EBB8F', 
           x = 455, y = 0.992, 
           label = "p=0.455", 
           hjust = 0, vjust = 0, size = 4.5) + 
  annotate("text", color = '#886A9E', 
           x = 455, y = 0.978, 
           label = "p<2e-16*", 
           hjust = 0, vjust = 0, size = 4.5)

# ggsave('./plot/compare_lineage_gam_model/gam_pred~distance&ANI_col2.pdf', width = 4, height = 3.5)



# ------------------------------
# Lineage 1.1 time_interval
# ------------------------------

dist_seq <- seq(
  min(fit_data_11$time_interval, na.rm = TRUE),
  max(fit_data_11$time_interval, na.rm = TRUE),
  length.out = 200
)

other_means <- data.frame(
  distance_km = mean(fit_data_11$distance_km, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_11$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_11$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_11$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_11$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

newdata_11 <- cbind(
  data.frame(time_interval = dist_seq),
  other_means
)

pred <- predict(gam_model_11_ani, newdata_11, type = "response", se.fit = TRUE)
newdata_11$pred  <- pred$fit
newdata_11$upper <- pred$fit + 1.96 * pred$se.fit
newdata_11$lower <- pred$fit - 1.96 * pred$se.fit

scale_val  <- scaling_params[scaling_params$variable == 'time_interval', 'scale']
center_val <- scaling_params[scaling_params$variable == 'time_interval', 'center']
newdata_11$time_interval <- newdata_11$time_interval * scale_val + center_val

# ------------------------------
# Lineage 1.3 time_interval
# ------------------------------

dist_seq <- seq(
  min(fit_data_13$time_interval, na.rm = TRUE),
  max(fit_data_13$time_interval, na.rm = TRUE),
  length.out = 200
)

other_means <- data.frame(
  distance_km = mean(fit_data_13$distance_km, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_13$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_13$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_13$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_13$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

newdata_13 <- cbind(
  data.frame(time_interval = dist_seq),
  other_means
)

pred <- predict(gam_model_13_ani, newdata_13, type = "response", se.fit = TRUE)
newdata_13$pred  <- pred$fit
newdata_13$upper <- pred$fit + 1.96 * pred$se.fit
newdata_13$lower <- pred$fit - 1.96 * pred$se.fit

scale_val  <- scaling_params[scaling_params$variable == 'time_interval', 'scale']
center_val <- scaling_params[scaling_params$variable == 'time_interval', 'center']
newdata_13$time_interval <- newdata_13$time_interval * scale_val + center_val



# ------------------------------
# Merge Plot time_interval
# ------------------------------

# p-value
summary(gam_model_11_ani) # s(time_interval) 0.253
summary(gam_model_13_ani) # s(time_interval) 0.371

p2 <- ggplot(data = NULL, aes(x = time_interval)) + 
  # geom_point(data = fit_data_11, aes(x = distance_km, y = ANI), alpha = 0.5) +
  geom_ribbon(data = newdata_11, aes(ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) +
  geom_line(data = newdata_11, aes(y = pred), color = '#7EBB8F', linewidth = 1) + 
  geom_line(data = newdata_11, aes(y = upper), color = '#7EBB8F', linetype = 2) +
  geom_line(data = newdata_11, aes(y = lower), color = '#7EBB8F', linetype = 2) + 
  geom_ribbon(data = newdata_13, aes(ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) + 
  geom_line(data = newdata_13, aes(y = pred), color = '#886A9E', linewidth = 1) +
  geom_line(data = newdata_13, aes(y = upper), color = '#886A9E', linetype = 2) +
  geom_line(data = newdata_13, aes(y = lower), color = '#886A9E', linetype = 2) +
  coord_cartesian(xlim = c(100, NA), ylim = c(0.975, 0.999)) + 
  theme_classic(base_size = 18) + 
  labs(
    x = "Time interval (days)",
    y = "ANI"
  ) + 
  annotate("text", color = '#7EBB8F', 
           x = 1000, y = 0.997, 
           label = "p=0.253", 
           hjust = 0, vjust = 0, size = 4.5) + 
  annotate("text", color = '#886A9E', 
           x = 1000, y = 0.985, 
           label = "p=0.371", 
           hjust = 0, vjust = 0, size = 4.5)

# ggsave('./plot/compare_lineage_gam_model/gam_pred~time&ANI_col2.pdf', width = 4, height = 3.5)



# ------------------------------
# Gene content similarity Lineage 1.1 distance_km
# ------------------------------

# 1. 生成 distance_km 等间隔序列
dist_seq <- seq(
  min(fit_data_11$distance_km, na.rm = TRUE),
  max(fit_data_11$distance_km, na.rm = TRUE),
  length.out = 200
)

# 2. 固定其他变量
#   如果模型中有其他连续变量，可设它们 = mean(...) 或 median(...)
#   如果有因子变量，可以选一个参照水平
other_means <- data.frame(
  time_interval = mean(fit_data_11$time_interval, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_11$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_11$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_11$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_11$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

# 3. 构造 newdata
newdata_11 <- cbind(
  data.frame(distance_km = dist_seq),
  other_means
)

# 4. 预测到原始尺度，并获取标准误
pred <- predict(gam_model_11_gc, newdata_11, type = "response", se.fit = TRUE)
newdata_11$pred  <- pred$fit
newdata_11$upper <- pred$fit + 1.96 * pred$se.fit
newdata_11$lower <- pred$fit - 1.96 * pred$se.fit

# 5. 将 distance_km 转换回原始尺度
scale_val  <- scaling_params[scaling_params$variable == 'distance_km', 'scale']
center_val <- scaling_params[scaling_params$variable == 'distance_km', 'center']
newdata_11$distance_km <- newdata_11$distance_km * scale_val + center_val

# ------------------------------
# Gene content similarity Lineage 1.3 distance_km
# ------------------------------

dist_seq <- seq(
  min(fit_data_13$distance_km, na.rm = TRUE),
  max(fit_data_13$distance_km, na.rm = TRUE),
  length.out = 200
)

other_means <- data.frame(
  time_interval = mean(fit_data_13$time_interval, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_13$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_13$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_13$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_13$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

newdata_13 <- cbind(
  data.frame(distance_km = dist_seq),
  other_means
)

pred <- predict(gam_model_13_gc, newdata_13, type = "response", se.fit = TRUE)
newdata_13$pred  <- pred$fit
newdata_13$upper <- pred$fit + 1.96 * pred$se.fit
newdata_13$lower <- pred$fit - 1.96 * pred$se.fit

scale_val  <- scaling_params[scaling_params$variable == 'distance_km', 'scale']
center_val <- scaling_params[scaling_params$variable == 'distance_km', 'center']
newdata_13$distance_km <- newdata_13$distance_km * scale_val + center_val



# ------------------------------
# Merge Plot distance_km
# ------------------------------

# p-value
summary(gam_model_11_gc) # s(distance_km) 0.501
summary(gam_model_13_gc) # s(distance_km) <2e-16

p3 <- ggplot() + 
  # geom_point(data = fit_data_11, aes(x = distance_km, y = ANI), alpha = 0.5) +
  geom_ribbon(data = newdata_11, aes(x = distance_km, ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) +
  geom_line(data = newdata_11, aes(x = distance_km, y = pred), color = '#7EBB8F', linewidth = 1) + 
  geom_line(data = newdata_11, aes(x = distance_km, y = upper), color = '#7EBB8F', linetype = 2) +
  geom_line(data = newdata_11, aes(x = distance_km, y = lower), color = '#7EBB8F', linetype = 2) + 
  geom_ribbon(data = newdata_13, aes(x = distance_km, ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) + 
  geom_line(data = newdata_13, aes(x = distance_km, y = pred), color = '#886A9E', linewidth = 1) +
  geom_line(data = newdata_13, aes(x = distance_km, y = upper), color = '#886A9E', linetype = 2) +
  geom_line(data = newdata_13, aes(x = distance_km, y = lower), color = '#886A9E', linetype = 2) +
  coord_cartesian(xlim = c(100, NA), ylim = c(0.89, 0.98)) +
  theme_classic(base_size = 18) + 
  labs(
    x = "Geographic distance (km)",
    y = "Gene content similarity"
  ) + 
  annotate("text", color = '#7EBB8F', 
           x = 460, y = 0.955, 
           label = "p=0.501", 
           hjust = 0, vjust = 0, size = 4.5) + 
  annotate("text", color = '#886A9E', 
           x = 460, y = 0.910, 
           label = "p<2e-16*", 
           hjust = 0, vjust = 0, size = 4.5)

# ggsave('./plot/compare_lineage_gam_model/gam_pred~distance&gene_content_col2.pdf', width = 4, height = 3.5)



# ------------------------------
# Gene content similarity Lineage 1.1 time_interval
# ------------------------------

dist_seq <- seq(
  min(fit_data_11$time_interval, na.rm = TRUE),
  max(fit_data_11$time_interval, na.rm = TRUE),
  length.out = 200
)

other_means <- data.frame(
  distance_km = mean(fit_data_11$distance_km, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_11$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_11$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_11$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_11$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

newdata_11 <- cbind(
  data.frame(time_interval = dist_seq),
  other_means
)

pred <- predict(gam_model_11_gc, newdata_11, type = "response", se.fit = TRUE)
newdata_11$pred  <- pred$fit
newdata_11$upper <- pred$fit + 1.96 * pred$se.fit
newdata_11$lower <- pred$fit - 1.96 * pred$se.fit

scale_val  <- scaling_params[scaling_params$variable == 'time_interval', 'scale']
center_val <- scaling_params[scaling_params$variable == 'time_interval', 'center']
newdata_11$time_interval <- newdata_11$time_interval * scale_val + center_val

# ------------------------------
# Gene content similarity Lineage 1.3 time_interval
# ------------------------------

dist_seq <- seq(
  min(fit_data_13$time_interval, na.rm = TRUE),
  max(fit_data_13$time_interval, na.rm = TRUE),
  length.out = 200
)

other_means <- data.frame(
  distance_km = mean(fit_data_13$distance_km, na.rm = TRUE),
  yearly_surface_pressure_2m_difference     = mean(fit_data_13$yearly_surface_pressure_2m_difference, na.rm = TRUE),
  yearly_dewpoint_temperature_2m_difference = mean(fit_data_13$yearly_dewpoint_temperature_2m_difference, na.rm = TRUE),
  yearly_temperature_2m_difference          = mean(fit_data_13$yearly_temperature_2m_difference, na.rm = TRUE),
  yearly_total_precipitation_sum_difference = mean(fit_data_13$yearly_total_precipitation_sum_difference, na.rm = TRUE)
)

newdata_13 <- cbind(
  data.frame(time_interval = dist_seq),
  other_means
)

pred <- predict(gam_model_13_gc, newdata_13, type = "response", se.fit = TRUE)
newdata_13$pred  <- pred$fit
newdata_13$upper <- pred$fit + 1.96 * pred$se.fit
newdata_13$lower <- pred$fit - 1.96 * pred$se.fit

scale_val  <- scaling_params[scaling_params$variable == 'time_interval', 'scale']
center_val <- scaling_params[scaling_params$variable == 'time_interval', 'center']
newdata_13$time_interval <- newdata_13$time_interval * scale_val + center_val



# ------------------------------
# Merge Plot time_interval
# ------------------------------

# p-value
summary(gam_model_11_gc) # s(time_interval) 0.630
summary(gam_model_13_gc) # s(time_interval) 0.476

p4 <- ggplot(data = NULL, aes(x = time_interval)) + 
  # geom_point(data = fit_data_11, aes(x = distance_km, y = ANI), alpha = 0.5) +
  geom_ribbon(data = newdata_11, aes(ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) +
  geom_line(data = newdata_11, aes(y = pred), color = '#7EBB8F', linewidth = 1) + 
  geom_line(data = newdata_11, aes(y = upper), color = '#7EBB8F', linetype = 2) +
  geom_line(data = newdata_11, aes(y = lower), color = '#7EBB8F', linetype = 2) + 
  geom_ribbon(data = newdata_13, aes(ymin = lower, ymax = upper), fill = 'grey75', alpha = 0.15) + 
  geom_line(data = newdata_13, aes(y = pred), color = '#886A9E', linewidth = 1) +
  geom_line(data = newdata_13, aes(y = upper), color = '#886A9E', linetype = 2) +
  geom_line(data = newdata_13, aes(y = lower), color = '#886A9E', linetype = 2) +
  coord_cartesian(xlim = c(100, NA), ylim = c(0.89, 0.98)) + 
  theme_classic(base_size = 18) + 
  labs(
    x = "Time interval (days)",
    y = "Gene content similarity"
  ) + 
  annotate("text", color = '#7EBB8F',
           x = 1000, y = 0.975,
           label = "p=0.253",
           hjust = 0, vjust = 0, size = 4.5) +
  annotate("text", color = '#886A9E',
           x = 1000, y = 0.925,
           label = "p=0.371",
           hjust = 0, vjust = 0, size = 4.5)

# ggsave('./plot/compare_lineage_gam_model/gam_pred~time&gene_content_col2.pdf', width = 4, height = 3.5)


cowplot::plot_grid(
  p1, p2, p3, p4, 
  ncol = 4, align = "v"
)

# ggsave('./plot/compare_lineage_gam_model/gam_pred_col2.pdf', width = 14, height = 3.5)



# p_obj <- plot(gam_model_11_ani, residuals = TRUE)
# p_obj <- p_obj[[1]] # distance_km
# sm_df   <- as.data.frame(p_obj[c("x", "se", "fit")])
# data_df <- as.data.frame(p_obj[c("raw", "p.resid")])
# 
# ggplot(sm_df, aes(x = x, y = fit)) +
#   geom_rug(data = data_df, mapping = aes(x = raw, y = NULL),
#            sides = "b") +
#   geom_point(data = data_df, mapping = aes(x = raw, y = p.resid)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
#               alpha = 0.3) +
#   geom_line() +
#   labs(x = p_obj$xlab, y = p_obj$ylab)


summary(gam_model_11_gc)
summary(gam_model_13_gc)

# 提取函数：包含 Estimate、StdError、p-value
extract_parametric <- function(model, model_name) {
  s <- summary(model)
  coefs <- as.data.frame(s$p.table)
  coefs$Variable <- rownames(coefs)
  coefs <- coefs %>%
    filter(Variable %in% c("yearly_surface_pressure_2m_difference",
                           "yearly_dewpoint_temperature_2m_difference",
                           "yearly_temperature_2m_difference",
                           "yearly_total_precipitation_sum_difference")) %>%
    mutate(
      Estimate = Estimate,
      StdError = `Std. Error`,
      p_value = `Pr(>|t|)`,
      Model = model_name,
      Lower = Estimate - 1.96 * StdError,
      Upper = Estimate + 1.96 * StdError,
      Significance = ifelse(p_value < 0.05, "Significant", "Not Significant")
    )
  return(coefs[, c("Variable", "Estimate", "StdError", "Lower", "Upper", "p_value", "Model", "Significance")])
}

# 提取两个模型的结果
df_model11 <- extract_parametric(gam_model_11_gc, "Lineage 1.1")
df_model13 <- extract_parametric(gam_model_13_gc, "Lineage 1.3")

# 合并数据
all_data <- bind_rows(df_model11, df_model13) %>%
  mutate(
    color_group = case_when(
      Significance == "Significant" & Model == "Lineage 1.1" ~ "Lineage 1.1",
      Significance == "Significant" & Model == "Lineage 1.3" ~ "Lineage 1.3",
      TRUE ~ "Not Significant"
    ),
    # 单侧误差棒
    errorbar_min = ifelse(Estimate < 0, Lower, 1e-5),
    errorbar_max = ifelse(Estimate < 0, -1e-5, Upper)
  )

ggplot(all_data, aes(x = Estimate, y = Variable, fill = color_group)) + 
  geom_errorbarh(aes(xmin = errorbar_min, xmax = errorbar_max, color = color_group),
                height = 0.2, linewidth = 1,
                position = position_dodge(width = 0.5)) + 
  geom_col(position = position_dodge(width = 0.5), width = 0.6, linewidth = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) + 
  labs(x = "Estimate (with one-sided 95% CI)", y = "Variable") +
  theme_bw(base_size = 18) + 
  scale_color_manual(
    values = c("Lineage 1.1" = "#151d29", 
               "Lineage 1.3" = "#151d29", 
               "Not Significant" = "grey85")
  ) + 
  scale_fill_manual(
    values = c("Lineage 1.1" = "#7EBB8F", 
               "Lineage 1.3" = "#886A9E", 
               "Not Significant" = "grey85"),
    name = "Model"
  ) +
  theme(legend.position = "top") +
  coord_flip(xlim = c(-0.005, 0.005), ylim = c(0.5, 4.5), expand = FALSE) + 
  facet_grid(~ Model) +
  # 自定义背景填色：手动给分面上色
  ggh4x::facet_wrap2(~ Model, strip = ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = c("Lineage 1.1" = "#7EBB8F", "Lineage 1.3" = "#886A9E"))
  )) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave('./plot/compare_lineage_gam_model/forest~gene_content_col2_type2.pdf', width = 4.5, height = 8)



# 提取两个模型的结果
df_model11 <- extract_parametric(gam_model_11_ani, "Lineage 1.1")
df_model13 <- extract_parametric(gam_model_13_ani, "Lineage 1.3")

# 合并数据
all_data <- bind_rows(df_model11, df_model13) %>%
  mutate(
    color_group = case_when(
      Significance == "Significant" & Model == "Lineage 1.1" ~ "Lineage 1.1",
      Significance == "Significant" & Model == "Lineage 1.3" ~ "Lineage 1.3",
      TRUE ~ "Not Significant"
    ),
    # 单侧误差棒
    errorbar_min = ifelse(Estimate < 0, Lower, 1e-5),
    errorbar_max = ifelse(Estimate < 0, -1e-5, Upper)
  )

ggplot(all_data, aes(x = Estimate, y = Variable, fill = color_group)) + 
  geom_errorbarh(aes(xmin = errorbar_min, xmax = errorbar_max, color = color_group),
                 height = 0.2, linewidth = 1,
                 position = position_dodge(width = 0.5)) + 
  geom_col(position = position_dodge(width = 0.5), width = 0.6, linewidth = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) + 
  labs(x = "Estimate (with one-sided 95% CI)", y = "Variable") +
  theme_bw(base_size = 18) + 
  scale_color_manual(
    values = c("Lineage 1.1" = "#151d29", 
               "Lineage 1.3" = "#151d29", 
               "Not Significant" = "grey85")
  ) + 
  scale_fill_manual(
    values = c("Lineage 1.1" = "#7EBB8F", 
               "Lineage 1.3" = "#886A9E", 
               "Not Significant" = "grey85"),
    name = "Model"
  ) +
  theme(legend.position = "top") +
  coord_flip(xlim = c(-0.0015, 0.0015), ylim = c(0.5, 4.5), expand = FALSE) + 
  facet_grid(~ Model) +
  # 自定义背景填色：手动给分面上色
  ggh4x::facet_wrap2(~ Model, strip = ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = c("Lineage 1.1" = "#7EBB8F", "Lineage 1.3" = "#886A9E"))
  )) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave('./plot/compare_lineage_gam_model/forest~ani_col2_type2.pdf', width = 4.5, height = 8)
