library(tidyverse)
library(Biostrings)
library(ape)
library(phangorn)
library(ggtree)
library(ggmsa)
library(ecodist)
library(future.apply)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

# load('./r_image/compare_lineage_yunnan_v2.RData')
load('./r_image/metadata_v3_with_pairwise.RData')

# 读取 codon 距离矩阵 list（每个基因包含 $p 和 $ani 矩阵）
codon_distances <- readRDS('./r_image/lineage_core_genes_codon_distances_p_and_ani_panaroo.RDS')

# predictor 自变量名
# predictor_cols <- c(
#   "distance_km",
#   "time_interval",
#   "yearly_surface_pressure_2m_difference",
#   "yearly_dewpoint_temperature_2m_difference",
#   "yearly_temperature_2m_difference",
#   "yearly_total_precipitation_sum_difference"
# )

predictor_cols <- c(
  "distance_km"
)

# 构造 predictor 的距离矩阵（以 tip_1、tip_2 为行列名）
# -> for Roary
# pivot_to_dist <- function(df, col) {
#   tip_names <- unique(c(df$tip_1, df$tip_2))
#   dist_mat <- matrix(NA, nrow = length(tip_names), ncol = length(tip_names),
#                      dimnames = list(tip_names, tip_names))
#   for (i in seq_len(nrow(df))) {
#     x <- df$tip_1[i]
#     y <- df$tip_2[i]
#     val <- df[[col]][i]
#     if (!is.na(val)) {
#       dist_mat[x, y] <- val
#       dist_mat[y, x] <- val
#     }
#   }
#   return(as.dist(dist_mat))
# }

# -> for Panaroo
pivot_to_dist <- function(df, col) {
  tip_names <- unique(c(df$tip_i, df$tip_j))
  dist_mat <- matrix(NA, nrow = length(tip_names), ncol = length(tip_names),
                     dimnames = list(tip_names, tip_names))
  for (i in seq_len(nrow(df))) {
    x <- df$tip_i[i]
    y <- df$tip_j[i]
    val <- df[[col]][i]
    if (!is.na(val)) {
      dist_mat[x, y] <- val
      dist_mat[y, x] <- val
    }
  }
  return(as.dist(dist_mat))
}

# 运行 MRM 分析
run_mrm_analysis <- function(response_mat,
                             fit_data,
                             predictor_cols,
                             nperm = 999) {
  response_dist <- as.dist(response_mat)
  dist_list <- lapply(predictor_cols, function(col) {
    pivot_to_dist(fit_data, col)
  })
  names(dist_list) <- predictor_cols
  
  pred_names_dist <- paste0(predictor_cols, "_dist")
  rhs <- paste(pred_names_dist, collapse = " + ")
  formula_str <- paste0("response_dist ~ ", rhs)
  
  env_list <- list(response_dist = response_dist)
  for (i in seq_along(predictor_cols)) {
    env_list[[ pred_names_dist[i] ]] <- dist_list[[i]]
  }
  tmp_env <- new.env()
  list2env(env_list, envir = tmp_env)
  
  formula_obj <- as.formula(formula_str, env = tmp_env)
  
  mrm_result <- with(
    data = tmp_env,
    MRM(
      formula     = formula_obj,
      data        = NULL,
      nperm       = nperm,
      mrank       = TRUE
    )
  )
  
  mrm_result
}

# 准备 pairwise 数据
# -> for Roary
# pairwise_df <- pairwise_sequenced_cases_data %>% 
#   dplyr::select(tip_pair, tip_1, tip_2, all_of(predictor_cols))

# -> for Panaroo
pairwise_df <- metadata.pairwise %>%
  dplyr::select(tip_pair, tip_i, tip_j, all_of(predictor_cols))



# 取消 future.globals.maxSize 默认的 500 MB 限制
options(future.globals.maxSize = +Inf)  # 或设置成 1GB，如 1024 * 1024^2

# 设置并行
plan(multisession)

# -------------------------------
# 1）p-distance 建模
# -------------------------------
p_mrm_results <- future_lapply(names(codon_distances), function(gene) {
  mat <- codon_distances[[gene]]$p
  if (is.null(mat)) return(NULL)
  
  mrm_out <- try(run_mrm_analysis(
    response_mat = mat,
    fit_data = pairwise_df,
    predictor_cols = predictor_cols,
    nperm = 999
  ), silent = TRUE)
  
  if (inherits(mrm_out, "try-error")) return(NULL)
  
  coef_table <- mrm_out$coef
  coef_table <- coef_table[rownames(coef_table) != "Int", , drop = FALSE]
  rownames(coef_table) <- gsub("_dist$", "", rownames(coef_table))
  
  coefs <- pvals <- rep(NA_real_, length(predictor_cols))
  names(coefs) <- names(pvals) <- predictor_cols
  
  for (var in predictor_cols) {
    if (var %in% rownames(coef_table)) {
      coefs[var] <- coef_table[var, "response_dist"]
      pvals[var] <- coef_table[var, "pval"]
    }
  }
  
  return(list(coef = coefs, pval = pvals))
})
names(p_mrm_results) <- names(codon_distances)

# saveRDS(p_mrm_results, './r_image/lineage_core_genes_p_mrm_results_panaroo.RDS')

# -------------------------------
# 2）1-ANI 建模
# -------------------------------
ani_mrm_results <- future_lapply(names(codon_distances), function(gene) {
  mat <- codon_distances[[gene]]$ani
  if (is.null(mat)) return(NULL)
  
  mrm_out <- try(run_mrm_analysis(
    response_mat = mat,
    fit_data = pairwise_df,
    predictor_cols = predictor_cols,
    nperm = 999
  ), silent = TRUE)
  
  if (inherits(mrm_out, "try-error")) return(NULL)
  
  coef_table <- mrm_out$coef
  coef_table <- coef_table[rownames(coef_table) != "Int", , drop = FALSE]
  rownames(coef_table) <- gsub("_dist$", "", rownames(coef_table))
  
  coefs <- pvals <- rep(NA_real_, length(predictor_cols))
  names(coefs) <- names(pvals) <- predictor_cols
  
  for (var in predictor_cols) {
    if (var %in% rownames(coef_table)) {
      coefs[var] <- coef_table[var, "response_dist"]
      pvals[var] <- coef_table[var, "pval"]
    }
  }
  
  return(list(coef = coefs, pval = pvals))
})
names(ani_mrm_results) <- names(codon_distances)

# saveRDS(ani_mrm_results, './r_image/lineage_core_genes_ani_mrm_results_panaroo.RDS')
