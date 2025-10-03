library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/compare_lineage_yunnan_v2.RData')



datmat <- read.table('./data/roary_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat <- t(datmat) %>% as.data.frame()

colnames(datmat) <- datmat[1, ]
datmat <- datmat[-1, ]

rownames_datamat <- datmat$Gene
datmat$Gene <- NULL
datmat <- datmat %>% reframe(across(everything(), as.integer))
rownames(datmat) <- rownames_datamat
rm(rownames_datamat)

metadata_tree <- metadata_tree %>% dplyr::filter(tip != "Reference")



calculate_pangenome_curves <- function(datmat, break_step = 5, rep = 100) {
  # 假定 datmat 为样本 x 基因的 0/1 矩阵，需转置为基因 x 样本
  samples_matrix <- t(datmat)
  # 保证数值类型（若 datmat 已经转换，可省略此步）
  samples_matrix <- as.data.frame(lapply(as.data.frame(samples_matrix), as.integer),
                                  row.names = rownames(samples_matrix))
  
  # 内置函数：根据起始值、最大值和步长确定采样断点
  lastnum_seq <- function(w, x, y) {
    seq_breaks <- seq(w, x, by = y)
    return(c(seq_breaks, x))
  }
  
  # 对给定样本数采样并计算核心与全基因组基因数
  sampling_and_counting_genes <- function(matrix, number) {
    # 随机抽取 number 个样本（列）
    subset <- matrix[, sample(ncol(matrix), size = number)]
    # 每个基因在选定样本中的出现次数
    subset <- cbind(subset, Total = rowSums(subset))
    # 核心基因：在所有样本中都有出现
    val_core <- sum(subset$Total == number)
    # 全基因组：在至少一个样本中出现
    val_pan <- sum(subset$Total > 0)
    vals <- c(val_core, val_pan, number)
    return(vals)
  }
  
  # 对某一采样数重复 rep 次，计算核心基因数的五数概括
  core_genome_values <- function(matrix, number, rep = rep) {
    rep_sample <- replicate(rep, sampling_and_counting_genes(matrix, number))
    core_stats <- quantile(rep_sample[1,])
    samples <- unique(rep_sample[3,])
    core_data <- as.data.frame(t(c(core_stats, "Number of samples" = samples)))
    return(core_data)
  }
  
  # 同理，计算全基因组（pan genome）基因数的五数概括
  pan_genome_values <- function(matrix, number, rep = rep) {
    rep_sample <- replicate(rep, sampling_and_counting_genes(matrix, number))
    pan_stats <- quantile(rep_sample[2,])
    samples <- unique(rep_sample[3,])
    pan_data <- as.data.frame(t(c(pan_stats, "Samples" = samples)))
    return(pan_data)
  }
  
  # 根据样本总数确定采样断点
  sample_num <- ncol(samples_matrix)
  sizes <- lastnum_seq(break_step, sample_num, break_step)
  
  # 初始化存储结果的数据框
  final_core_curve_df <- data.frame()
  final_pan_curve_df <- data.frame()
  
  # 对每个采样断点进行 rep 次采样并存储统计量
  for (n in sizes) {
    core_values <- core_genome_values(samples_matrix, n, rep = rep)
    final_core_curve_df <- rbind(final_core_curve_df, core_values)
    
    pan_values <- pan_genome_values(samples_matrix, n, rep = rep)
    final_pan_curve_df <- rbind(final_pan_curve_df, pan_values)
  }
  
  # 整理列名
  colnames(final_core_curve_df) <- c("Minimum", "First_Quantile", "Median", "Third_Quantile", "Maximum", "Samples")
  colnames(final_pan_curve_df)  <- c("Minimum", "First_Quantile", "Median", "Third_Quantile", "Maximum", "Samples")
  
  # 返回CORE_CURVE和PAN_CURVE
  return(list(CORE_CURVE = final_core_curve_df, PAN_CURVE = final_pan_curve_df))
}

plotNewGenesReg <- function(pan, col = "firebrick") {
  
  ggplotRegression <- function (fit,my_col,xlb,ylb) {
    require(ggplot2)
    # 提取模型中自变量和因变量的名称
    xvar <- names(fit$model)[2]
    yvar <- names(fit$model)[1]
    ggplot(fit$model, aes(x = !!sym(xvar), y = !!sym(yvar))) + 
      geom_point(col = my_col) +
      stat_smooth(method = "lm", col = my_col) + 
      xlab(xlb) + ylab(ylb) + 
      labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                         " Slope =",signif(fit$coef[[2]], 5),
                         "Intercept =",signif(fit$coef[[1]],5 ),
                         " P =",signif(summary(fit)$coef[2,4], 5))) + 
      theme_bw(base_size = 14)
  }
  
  ##### creo dataframe con le differenze tra le mediane al crescere dei genomi e il numero di genomi samplati
  m <- pan$Median
  d <- m[2:length(m)] - m[1:length(m) -1]
  df <- data.frame(d,pan$Samples[2:length(pan$Samples)])
  colnames(df) <- c("New_genes","N")
  
  ##### uso funzione che mi crea plot e aggiunge stats di lm sopra
  dfl <- log10(df)
  
  ### remove any NA / -Inf due to log 
  dfl <- dfl[!is.na(dfl$New_genes),]
  dfl <- dfl[dfl$New_genes != "-Inf",]
  
  ### plot regression
  g <- ggplotRegression(lm(New_genes ~ N, data = dfl), col,"Log Number of genomes", "Log Number of new genes")
  
  g
}


#######################################
## Lineage 1.1 、1.3 的基因开放性检测


# view(curves_11$PAN_CURVE)
# view(curves_13$CORE_CURVE)

length(metadata_tree$tip[metadata_tree$group_1 == '1'])
# [1] 657
curves_L1 <- calculate_pangenome_curves(datmat[metadata_tree$tip[metadata_tree$group_1 == '1'], ], break_step = 10, rep = 1000)
plotNewGenesReg(curves_L1$PAN_CURVE)

length(metadata_tree$tip[metadata_tree$group_1 == '2'])
# [1] 48
curves_L2 <- calculate_pangenome_curves(datmat[metadata_tree$tip[metadata_tree$group_1 == '2'], ], break_step = 5, rep = 1000)
plotNewGenesReg(curves_L2$PAN_CURVE)

length(metadata_tree$tip[metadata_tree$group_1 == '3'])
# [1] 140
curves_L3 <- calculate_pangenome_curves(datmat[metadata_tree$tip[metadata_tree$group_1 == '3'], ], break_step = 5, rep = 1000)
plotNewGenesReg(curves_L3$PAN_CURVE)

# ggsave('./plot/pan_genome_analysis/open_close_genome_all.pdf', width = 4.5, height = 3.5)

curves_L1_grop <- curves_L1
curves_L1_grop$PAN_CURVE <- curves_L1_grop$PAN_CURVE[1:(nrow(curves_L1_grop$PAN_CURVE) - 1), ]
curves_L2_grop <- curves_L2
curves_L2_grop$PAN_CURVE <- curves_L2_grop$PAN_CURVE[1:(nrow(curves_L2_grop$PAN_CURVE) - 1), ]
curves_L3_grop <- curves_L3
curves_L3_grop$PAN_CURVE <- curves_L3_grop$PAN_CURVE[1:(nrow(curves_L3_grop$PAN_CURVE) - 1), ]

# 计算新基因数：相邻中位数之差（去除第一个点），并记录对应样本数
mL1 <- curves_L1_grop$PAN_CURVE$Median
NL1 <- curves_L1_grop$PAN_CURVE$Samples
dL1 <- mL1[-1] - mL1[-length(mL1)]
dfL1 <- data.frame(New_genes = dL1, N = NL1[-1], group = "1")

mL2 <- curves_L2_grop$PAN_CURVE$Median
NL2 <- curves_L2_grop$PAN_CURVE$Samples
dL2 <- mL2[-1] - mL2[-length(mL2)]
dfL2 <- data.frame(New_genes = dL2, N = NL2[-1], group = "2")

mL3 <- curves_L3_grop$PAN_CURVE$Median
NL3 <- curves_L3_grop$PAN_CURVE$Samples
dL3 <- mL3[-1] - mL3[-length(mL3)]
dfL3 <- data.frame(New_genes = dL3, N = NL3[-1], group = "3")

# 合并两组数据，并进行 log10 转换
df_combined <- rbind(dfL1, dfL2, dfL3)
df_log <- data.frame(
  New_genes = log10(df_combined$New_genes),
  N = log10(df_combined$N),
  group = df_combined$group
)
df_log <- df_log[is.finite(df_log$New_genes) & is.finite(df_log$N), ]

# 分别拟合两组数据的线性回归模型
fit_L1 <- lm(New_genes ~ N, data = subset(df_log, group == "1"))
fit_L2 <- lm(New_genes ~ N, data = subset(df_log, group == "2"))
fit_L3 <- lm(New_genes ~ N, data = subset(df_log, group == "3"))

# 构造显示回归统计信息的标签
label_L1 <- paste("Lineage 1:",
                  # "Adj R2 =", signif(summary(fit_L1)$adj.r.squared, 3),
                  # "Intercept =", signif(fit_L1$coef[[1]], 3),
                  "Slope =", signif(fit_L1$coef[[2]], 3),
                  "P =", signif(summary(fit_L1)$coef[2, 4], 3))
print(label_L1)
label_L2 <- paste("Lineage 2:",
                  # "Adj R2 =", signif(summary(fit_L2)$adj.r.squared, 3),
                  # "Intercept =", signif(fit_L2$coef[[1]], 3),
                  "Slope =", signif(fit_L2$coef[[2]], 3),
                  "P =", signif(summary(fit_L2)$coef[2, 4], 3))
print(label_L2)
label_L3 <- paste("Lineage 3:",
                  # "Adj R2 =", signif(summary(fit_L3)$adj.r.squared, 3),
                  # "Intercept =", signif(fit_L3$coef[[1]], 3),
                  "Slope =", signif(fit_L3$coef[[2]], 3),
                  "P =", signif(summary(fit_L3)$coef[2, 4], 3))
print(label_L3)

# 绘图：散点 + 回归线（含置信区间），并添加文本注释显示回归统计信息

# ggplot(df_log, aes(x = N, y = New_genes, color = group)) +
#   geom_point() +
#   stat_smooth(method = "lm", se = TRUE, fill = 'grey75', alpha = 0.15) +
#   scale_color_manual(values = c("1" = "#9DD0CF", "2" = "#D0A341", "3" = "#c9a7b7")) +
#   xlab("Log Number of genomes") +
#   ylab("Log Number of new genes") +
#   theme_bw(base_size = 14) + 
#   annotation_logticks() + 
#   coord_fixed(xlim = c(0.55, 2.1), ylim = c(0.55, 2.1))


# 定义一个函数用于生成预测数据，
# x轴从 1 到对应组数据中 N 的最大值，生成 n_points 个点
create_prediction_df <- function(model, group_data, group_label, n_points = 100) {
  x_min <- 1
  x_max <- max(group_data$N, na.rm = TRUE)
  newdata <- data.frame(N = seq(x_min, x_max, length.out = n_points))
  pred <- predict(model, newdata = newdata, se.fit = TRUE)
  newdata$fit <- pred$fit
  newdata$lwr <- pred$fit - 1.96 * pred$se.fit
  newdata$upr <- pred$fit + 1.96 * pred$se.fit
  newdata$group <- group_label
  return(newdata)
}

# 为每个组构建预测数据
pred1 <- create_prediction_df(fit_L1, subset(df_log, group == "1"), "1")
pred2 <- create_prediction_df(fit_L2, subset(df_log, group == "2"), "2")
pred3 <- create_prediction_df(fit_L3, subset(df_log, group == "3"), "3")
pred_df <- rbind(pred1, pred2, pred3)

# 选择合适的文本位置，根据数据范围稍作偏移
x_range <- range(df_log$N)
y_range <- range(df_log$New_genes)
x_pos <- x_range[1] + 0.2 * diff(x_range)
y_pos_L1 <- y_range[1] + 0.15 * diff(y_range)
y_pos_L2 <- y_pos_L1 - 0.1 * diff(y_range)
y_pos_L3 <- y_pos_L2 - 0.1 * diff(y_range)

ggplot() + 
  geom_point(data = df_log, aes(x = N, y = New_genes, color = group)) +
  geom_ribbon(data = pred_df, aes(x = N, ymin = lwr, ymax = upr, fill = group), alpha = 0.15, color = NA) +
  scale_color_manual(values = c("1" = "#9DD0CF", "2" = "#D0A341", "3" = "#c9a7b7")) +
  scale_fill_manual(values = c("1" = "#9DD0CF", "2" = "#D0A341", "3" = "#c9a7b7")) + 
  ggnewscale::new_scale_color() + 
  geom_line(data = pred_df, aes(x = N, y = fit, color = group), linewidth = 1) + 
  scale_color_manual(values = c("1" = "#38c2bf", "2" = "#c4901d", "3" = "#b86c90")) +
  xlab("Log Number of genomes") +
  ylab("Log Number of new genes") +
  theme_bw(base_size = 14) + 
  annotation_logticks() + 
  # coord_fixed(xlim = c(1, 3), ylim = c(1, 3)) 
  coord_cartesian(xlim = c(1, 3), ylim = c(1, 3.3)) + 
  annotate("text", x = x_pos, y = y_pos_L1, label = label_L1, 
           color = "black", hjust = 0, size = 3.5) +
  annotate("text", x = x_pos, y = y_pos_L2, label = label_L2, 
           color = "black", hjust = 0, size = 3.5) + 
  annotate("text", x = x_pos, y = y_pos_L3, label = label_L3, 
           color = "black", hjust = 0, size = 3.5)

# ggsave('./plot/pan_genome_analysis/open_close_genome_global_col2.pdf', width = 5, height = 4)

# 1. 确保group为因子
df_log$group <- as.factor(df_log$group)
# 2. 建立全模型：考虑组别和交互项（截距和斜率均可不同）
fit_full <- lm(New_genes ~ N * group, data = df_log)
# 3. 建立加性模型（只有截距可以不同，斜率相同）
fit_add <- lm(New_genes ~ N + group, data = df_log)
# 4. 模型比较：检验加入交互项是否显著改善模型拟合
anova_result <- anova(fit_add, fit_full)
anova_result

# Analysis of Variance Table
# 
# Model 1: New_genes ~ N + group
# Model 2: New_genes ~ N * group
#   Res.Df    RSS Df Sum of Sq      F Pr(>F)
# 1     62 1.2201                           
# 2     60 1.2164  2 0.0036515 0.0901  0.914


# save.image('./r_image/pan_genome_openess_global.RData')
