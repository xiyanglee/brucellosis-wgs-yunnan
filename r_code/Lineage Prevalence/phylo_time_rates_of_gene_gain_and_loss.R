# withr::with_libpaths("/public/r_share_library/4.1", remotes::install_github("gtonkinhill/panstripe", build_vignettes = TRUE))
library(panstripe)
library(ape)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')
source('./r_code/Lineage Prevalence/lib.R')



beast_tree <- treeio::read.beast('./data/beast_gubbins_contig/run_1/beast_GTRIG_clean_core_filtered.tree')

normalize_datamat <- function(datmat) {
  colnames(datmat) <- datmat[1, ]
  datmat <- datmat[-1, ]
  
  rownames_datamat <- datmat$Gene
  datmat$Gene <- NULL
  
  datmat <- datmat %>% reframe(across(everything(), as.integer))
  rownames(datmat) <- rownames_datamat
  datmat <- as.data.frame(t(datmat))
  
  return(datmat)
}

datmat.onehot <- read.table('./data/panaroo_result_sensitive/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat.onehot <- normalize_datamat(datmat.onehot)


get_anc_changes <- function(pa, tree) {
  # 对齐行顺序
  pa <- pa[match(tree$tip.label, rownames(pa)), , drop = FALSE]
  # 只保留在不同样本间有变化的基因
  index <- which(apply(pa, 2, function(x) length(unique(x))) > 1)
  
  # 一次性对所有基因调用 ASR
  res_list <- purrr::map(index, ~{
    asr_max_parsimony(tree, pa[, .x] + 1, Nstates = 2)
  })
  
  # 边的索引，用于抽取 changes_* 对应每条边
  edges <- tree$edge[,2]
  n_edges <- nrow(tree$edge)
  
  # 分别构建 gain 与 loss 矩阵：行=边，列=基因
  anc_gain <- vapply(res_list,
                     function(res) res$changes_gain[edges],
                     numeric(n_edges))
  anc_loss <- vapply(res_list,
                     function(res) res$changes_loss[edges],
                     numeric(n_edges))
  
  # 恢复列名，便于追踪原基因
  colnames(anc_gain) <- colnames(anc_loss) <- colnames(pa)[index]
  
  list(
    gain = anc_gain,   # 每条边上 0→1（获得）的期望次数/概率
    loss = anc_loss    # 每条边上 1→0（丢失）的期望次数/概率
  )
}

get_fit_data <- function(pa, tree, root_depth = 0) {
  anc_states <- get_anc_changes(pa, tree)
  
  dat <- tibble::tibble(
    acc_gain = rowSums(anc_states$gain, na.rm = TRUE),
    acc_loss = rowSums(anc_states$loss, na.rm = TRUE),
    core=tree$edge.length,
    istip=as.numeric(tree$edge[,2]<=length(tree$tip.label) )
  )
  
  # add depth and filter if requested
  dat$depth <- ape::node.depth.edgelength(tree)[tree$edge[,1]] + root_depth
  
  return(dat)
}

get_sub_fit_data <- function(pa, tree, tips_to_keep, use_root_depth = TRUE) {
  # Step 1: 找出子树根节点（即 tips_to_keep 的 MRCA）
  mrca_node <- getMRCA(tree, tips_to_keep)
  
  # Step 2: 计算原始树中该节点的 depth
  # 注意 ape::node.depth.edgelength 是从 tip 向 root 方向计算的，每个 node 到 tip 的总长
  # 但你要从 root 到 node 的路径，可用 ape::node.depth.edgelength(reverse = TRUE)
  node_depths <- node.depth.edgelength(tree)
  root_depth <- node_depths[mrca_node]
  
  # Step 3: drop.tip 提取子树
  subtree <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
  subpa   <- pa[subtree$tip.label, ]
  
  if(use_root_depth) {
    return(get_fit_data(subpa, subtree, root_depth = root_depth))
  } else {
    return(get_fit_data(subpa, subtree, root_depth = 0))
  }
}


## 全球样本 Lineage 1、2、3  的模型

# pa <- datmat.onehot
# tree <- beast_tree@phylo
# 
# # 要保留的tip标签
# tips_to_keep <- intersect(metadata.phylo$assembly_accession[metadata.phylo$lineage_level_1 == "1"], tree$tip.label)
# dat_l1 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1")
# 
# tips_to_keep <- intersect(metadata.phylo$assembly_accession[metadata.phylo$lineage_level_1 == "2"], tree$tip.label)
# dat_l2 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "2")
# 
# tips_to_keep <- intersect(metadata.phylo$assembly_accession[metadata.phylo$lineage_level_1 == "3"], tree$tip.label)
# dat_l3 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "3")
# 
# 
# dat <- rbind(dat_l1, dat_l2, dat_l3)
# # dat <- dplyr::filter(dat, depth > 2500)
# dat$lineage <- as.factor(dat$lineage)
# 
# ggplot(dat, aes(x = depth, y = acc_gain)) + 
#   geom_point(aes(color = lineage)) + 
#   theme_bw()
# 
# ggplot(dat, aes(x = depth, y = acc_loss)) + 
#   geom_point(aes(color = lineage)) + 
#   theme_bw()
# 
# m_gain <- stats::glm(acc_gain ~ 0 + core + depth + lineage:core, data = dat, family = "quasipoisson")
# m_loss <- stats::glm(acc_loss ~ 0 + core + depth + lineage:core, data = dat, family = "quasipoisson")
# 
# # m_gain <- fit_tweedie(acc_gain ~ 0 + depth + core + lineage:core, data = dat, method = 'glm')
# # m_loss <- fit_tweedie(acc_loss ~ 0 + depth + core + lineage:core, data = dat, method = 'glm')
# summary(m_gain)
# summary(m_loss)
# 
# marginaleffects::slopes(m_gain, variables = "core", by = "lineage") %>% 
#   as.data.frame() %>% 
#   ggplot(aes(x = estimate, y = lineage)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#   labs(
#     x = "margin effect",
#     y = "Lineage"
#   ) +
#   theme_minimal()
# 
# marginaleffects::slopes(m_loss, variables = "core", by = "lineage") %>% 
#   as.data.frame() %>% 
#   ggplot(aes(x = estimate, y = lineage)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#   labs(
#     x = "margin effect",
#     y = "Lineage"
#   ) +
#   theme_minimal()


## 云南样本 Lineage 1.1、1.3  模型

softgene_names <- colnames(datmat.onehot)[colMeans(datmat.onehot) > 0.15 & colMeans(datmat.onehot) < 0.95]
raregene_names <- colnames(datmat.onehot)[colMeans(datmat.onehot) <= 0.15]

################################################################################
## -> genome
pa <- datmat.onehot
tree <- beast_tree@phylo

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.1"], tree$tip.label)
dat_l11 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.1")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.3"], tree$tip.label)
dat_l13 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.3")


dat_genome <- rbind(dat_l11, dat_l13)
# dat_genome <- dplyr::filter(dat_genome, depth > 2500)
dat_genome$lineage <- as.factor(dat_genome$lineage)

ggplot(dat_genome, aes(x = log10(core))) + 
  geom_point(aes(y = acc_gain, color = lineage)) + 
  geom_point(aes(y = - acc_loss, color = lineage)) +
  theme_bw()

m_gain <- stats::glm(acc_gain ~ 0 + depth + core + lineage:core, data = dat_genome, family = "quasipoisson")
m_loss <- stats::glm(acc_loss ~ 0 + depth + core + lineage:core, data = dat_genome, family = "quasipoisson")
# m_gain <- fit_tweedie(acc_gain ~ 0 + depth + core + lineage:core + lineage:depth, data = dat_genome, method = 'glm')
# m_loss <- fit_tweedie(acc_loss ~ 0 + depth + core + lineage:core + lineage:depth, data = dat_genome, method = 'glm')
summary(m_gain)
# Call:
# stats::glm(formula = acc_gain ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_genome)
# 
# Deviance Residuals: 
#    Min      1Q  Median      3Q     Max  
# -6.344  -2.295  -0.899   1.246  15.994  
# 
# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0030547  0.0003874   7.885 2.45e-13 ***
# core            0.0016976  0.0033148   0.512  0.60916    
# core:lineage1.1 0.0122989  0.0039086   3.147  0.00192 ** 
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasipoisson family taken to be 46.07657)
# 
#     Null deviance: 5747.5  on 192  degrees of freedom
# Residual deviance: 3132.8  on 189  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 7
summary(m_loss)
# Call:
# stats::glm(formula = acc_loss ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_genome)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -6.2017  -2.2858  -0.8024   1.1209  15.4694  
# 
# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0029505  0.0004038   7.307 7.47e-12 ***
# core            0.0017793  0.0032489   0.548  0.58457    
# core:lineage1.1 0.0121472  0.0038724   3.137  0.00198 ** 
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasipoisson family taken to be 46.4699)
# 
#     Null deviance: 5478.0  on 192  degrees of freedom
# Residual deviance: 3085.6  on 189  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 7


# install.packages("marginaleffects", lib = "/public/r_share_library/4.1")
marginaleffects::plot_predictions(m_gain, condition = list("core", "lineage")) + 
  geom_point(data = dat_genome, aes(x = core, y = acc_gain, colour = lineage)) + 
  scale_y_log10() +
  # scale_x_log10() + 
  coord_cartesian(xlim = c(1, 300), ylim = c(1, 1000))

marginaleffects::slopes(m_gain, variables = "core", by = "lineage") %>% 
  as.data.frame() %>% 
  ggplot(aes(x = estimate, y = lineage)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "margin effect",
    y = "Lineage"
  ) +
  theme_minimal()

marginaleffects::slopes(m_loss, variables = "core", by = "lineage") %>% 
  as.data.frame() %>% 
  ggplot(aes(x = estimate, y = lineage)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "margin effect",
    y = "Lineage"
  ) +
  theme_minimal()

genome_margin_effect <- rbind(
  marginaleffects::slopes(m_gain, variables = "core", by = "lineage") %>% 
    as.data.frame() %>% 
    mutate(
      events = 'gain events', 
      gene_type = 'genome'
    ), 
  marginaleffects::slopes(m_loss, variables = "core", by = "lineage") %>% 
    as.data.frame() %>% 
    mutate(
      events = 'loss events', 
      gene_type = 'genome'
    )
)

genome_margin_effect
#   term contrast lineage   estimate  std.error statistic     p.value   s.value    conf.low  conf.high      events gene_type
# 1 core    dY/dX     1.1 0.05617317 0.01781030 3.1539711 0.001610651 9.2781407  0.02126563 0.09108072 gain events    genome
# 2 core    dY/dX     1.3 0.01991400 0.04035514 0.4934688 0.621681376 0.6857527 -0.05918062 0.09900862 gain events    genome
# 3 core    dY/dX     1.1 0.05425211 0.01757919 3.0861562 0.002027623 8.9459949  0.01979754 0.08870668 loss events    genome
# 4 core    dY/dX     1.3 0.01925823 0.03675608 0.5239467 0.600315589 0.7362070 -0.05278236 0.09129881 loss events    genome


################################################################################
## -> soft-gene
pa <- datmat.onehot[, softgene_names]
tree <- beast_tree@phylo

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.1"], tree$tip.label)
dat_l11 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.1")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.3"], tree$tip.label)
dat_l13 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.3")


dat_soft <- rbind(dat_l11, dat_l13)
dat_soft$lineage <- as.factor(dat_soft$lineage)

m_gain <- stats::glm(acc_gain ~ 0 + depth + core + lineage:core, data = dat_soft, family = "quasipoisson")
m_loss <- stats::glm(acc_loss ~ 0 + depth + core + lineage:core, data = dat_soft, family = "quasipoisson")
summary(m_gain)
# Call:
# stats::glm(formula = acc_gain ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_soft)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -4.7186  -2.2160  -0.9776   0.8358  15.4152  
# 
# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0027206  0.0003649   7.456 3.13e-12 ***
# core            0.0015455  0.0029292   0.528  0.59837    
# core:lineage1.1 0.0104054  0.0037087   2.806  0.00555 ** 
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasipoisson family taken to be 31.71713)
# 
#     Null deviance: 3933.1  on 192  degrees of freedom
# Residual deviance: 2372.3  on 189  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 7
summary(m_loss)
# Call:
# stats::glm(formula = acc_loss ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_soft)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -5.1264  -2.2099  -1.0760   0.5781  15.7923  
# 
# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0026782  0.0003973   6.740 1.86e-10 ***
# core            0.0016171  0.0030791   0.525  0.60007    
# core:lineage1.1 0.0118120  0.0036796   3.210  0.00156 ** 
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasipoisson family taken to be 36.72187)
# 
#     Null deviance: 4418.0  on 192  degrees of freedom
# Residual deviance: 2611.6  on 189  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 7

soft_gene_margin_effect <- rbind(
  marginaleffects::slopes(m_gain, variables = "core", by = "lineage") %>% 
    as.data.frame() %>% 
    mutate(
      events = 'gain events', 
      gene_type = 'soft_gene'
    ), 
  marginaleffects::slopes(m_loss, variables = "core", by = "lineage") %>% 
    as.data.frame() %>% 
    mutate(
      events = 'loss events', 
      gene_type = 'soft_gene'
    )
)

soft_gene_margin_effect
#   term contrast lineage   estimate  std.error statistic     p.value   s.value    conf.low  conf.high      events gene_type
# 1 core    dY/dX     1.1 0.03852339 0.01292583 2.9803412 0.002879275 8.4400789  0.01318923 0.06385756 gain events soft_gene
# 2 core    dY/dX     1.3 0.01383447 0.02722517 0.5081499 0.611348215 0.7099337 -0.03952588 0.06719482 gain events soft_gene
# 3 core    dY/dX     1.1 0.04726404 0.01468287 3.2189924 0.001286419 9.6024236  0.01848615 0.07604193 loss events soft_gene
# 4 core    dY/dX     1.3 0.01402713 0.02787568 0.5032031 0.614821504 0.7017605 -0.04060820 0.06866245 loss events soft_gene


################################################################################
## -> rare-gene
pa <- datmat.onehot[, raregene_names]
tree <- beast_tree@phylo

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.1"], tree$tip.label)
dat_l11 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.1")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.3"], tree$tip.label)
dat_l13 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.3")

dat_rare <- rbind(dat_l11, dat_l13)
dat_rare$lineage <- as.factor(dat_rare$lineage)

m_gain <- stats::glm(acc_gain ~ 0 + depth + core + lineage:core, data = dat_rare, family = "quasipoisson")
m_loss <- stats::glm(acc_loss ~ 0 + depth + core + lineage:core, data = dat_rare, family = "quasipoisson")

summary(m_gain)
# Call:
# stats::glm(formula = acc_gain ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_rare)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -2.8946  -1.4878  -1.4433  -0.0558  12.0498  
# 
# Coefficients: (1 not defined because of singularities)
#                   Estimate Std. Error t value Pr(>|t|)
# depth            0.0000884  0.0006185   0.143    0.886
# core            -0.0022660  0.0078102  -0.290    0.772
# core:lineage1.1  0.0129253  0.0080219   1.611    0.109
# core:lineage1.3         NA         NA      NA       NA
# 
# (Dispersion parameter for quasipoisson family taken to be 14.52837)
# 
#     Null deviance: 1198.52  on 192  degrees of freedom
# Residual deviance:  986.39  on 189  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 7
summary(m_loss)
# Call:
# stats::glm(formula = acc_loss ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_rare)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -1.9039  -1.4048  -1.3879  -0.5213  12.0270  
# 
# Coefficients: (1 not defined because of singularities)
#                   Estimate Std. Error t value Pr(>|t|)
# depth           -0.0001180  0.0006584  -0.179    0.858
# core            -0.0005967  0.0043895  -0.136    0.892
# core:lineage1.1  0.0039695  0.0064044   0.620    0.536
# core:lineage1.3         NA         NA      NA       NA
# 
# (Dispersion parameter for quasipoisson family taken to be 14.02131)
# 
#     Null deviance: 858.41  on 192  degrees of freedom
# Residual deviance: 852.21  on 189  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 8

rare_gene_margin_effect <- rbind(
  marginaleffects::slopes(m_gain, variables = "core", by = "lineage") %>% 
    as.data.frame() %>% 
    mutate(
      events = 'gain events', 
      gene_type = 'rare_gene'
    ), 
  marginaleffects::slopes(m_loss, variables = "core", by = "lineage") %>% 
    as.data.frame() %>% 
    mutate(
      events = 'loss events', 
      gene_type = 'rare_gene'
    )
)

rare_gene_margin_effect
#   term contrast lineage      estimate   std.error  statistic     p.value   s.value     conf.low   conf.high      events gene_type
# 1 core    dY/dX     1.1  0.0166457602 0.005686667  2.9271556 0.003420777 8.1914602  0.005500097 0.027791423 gain events rare_gene
# 2 core    dY/dX     1.3 -0.0022346756 0.007521823 -0.2970923 0.766396045 0.3838380 -0.016977177 0.012507826 gain events rare_gene
# 3 core    dY/dX     1.1  0.0035514320 0.005401229  0.6575230 0.510844701 0.9690433 -0.007034782 0.014137646 loss events rare_gene
# 4 core    dY/dX     1.3 -0.0005281402 0.003809003 -0.1386558 0.889722168 0.1685732 -0.007993648 0.006937368 loss events rare_gene



rbind(
  genome_margin_effect, 
  soft_gene_margin_effect, 
  rare_gene_margin_effect
) %>% 
  mutate(
    gene_type = factor(gene_type, levels = c("genome", "soft_gene", "rare_gene"))
  ) %>% 
  ggplot(aes(x = lineage, y = estimate, ymin = conf.low, ymax = conf.high, color = events)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +         # 中心点 + 置信区间线
  geom_hline(yintercept = 0, linetype = "dashed") +    # 零参考线
  scale_color_manual(values = c("gain events" = "#164a84", "loss events" = "#e2041b")) +
  labs(
    x = "Lineage",
    y = "Marginal Effect\n(β ± 95% CI)",
    color = "Event Type"
  ) + 
  facet_wrap(~ gene_type) + 
  theme_classic(base_size = 14)

# ggsave('./plot/Lineage Prevalence/phylo_time_rates_of_gene_gain_and_loss.pdf', width = 9, height = 3.1)

# save.image('./r_image/phylo_time_rates_of_gene_gain_and_loss.RData')


tmp_plot <- rbind(
  soft_gene_margin_effect, 
  rare_gene_margin_effect
) %>% 
  mutate(
    gene_type = factor(gene_type, levels = c("genome", "soft_gene", "rare_gene")), 
    sig_group = ifelse(p.value < 0.05, "p < 0.05", "p ≥ 0.05"),
    sig_group = factor(sig_group, levels = c("p < 0.05", "p ≥ 0.05"))
  )

p2 <- ggplot(tmp_plot, aes(x = gene_type, y = estimate, ymin = conf.low, ymax = conf.high, color = events)) + 
  geom_errorbar(position = position_dodge(width = 0.6), linewidth = 0.8, width = 0.2) + 
  geom_point(aes(shape = sig_group), position = position_dodge(width = 0.6), size = 3.5, stroke = 1, fill = "white") + 
  # geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("gain events" = "#164a84", "loss events" = "#e2041b")) + 
  scale_shape_manual(values = c("p ≥ 0.05" = 21, "p < 0.05" = 16)) +  # 16: 实心圆, 21: 空心圆
  labs(
    x = "Accessory gene type",
    y = "Rates of gene gain/loss\n(β ± 95% CI)",
    color = "Event Type"
  ) +
  facet_wrap(~ lineage) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),  # 去掉默认条
    strip.text = element_text(color = "white", face = "bold", size = 12), 
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  # 自定义背景填色：手动给分面上色
  ggh4x::facet_wrap2(~ lineage, strip = ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = c("1.1" = "#7EBB8F", "1.3" = "#886A9E"))
  )) +
  coord_cartesian(ylim = c(-0.055, 0.105)) + 
  theme(panel.background = element_rect(fill = "#f9f9f9"))

tmp_plot <- genome_margin_effect %>% 
  mutate(
    gene_type = factor(gene_type, levels = c("genome", "soft_gene", "rare_gene")), 
    sig_group = ifelse(p.value < 0.05, "p < 0.05", "p ≥ 0.05"),
    sig_group = factor(sig_group, levels = c("p < 0.05", "p ≥ 0.05"))
  )
  
p1 <- ggplot(tmp_plot, aes(x = lineage, y = estimate, ymin = conf.low, ymax = conf.high, color = events)) + 
  geom_errorbar(position = position_dodge(width = 0.6), linewidth = 0.8, width = 0.2) + 
  geom_point(aes(shape = sig_group), position = position_dodge(width = 0.6), size = 3.5, stroke = 1, fill = "white") + 
  # geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("gain events" = "#164a84", "loss events" = "#e2041b")) + 
  scale_shape_manual(values = c("p ≥ 0.05" = 21, "p < 0.05" = 16)) +  # 16: 实心圆, 21: 空心圆
  labs(
    x = "Lineage",
    y = "Rates of gene gain/loss\n(β ± 95% CI)",
    color = "Event Type"
  ) +
  facet_wrap(~ gene_type) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1),
    # strip.background = element_blank(),  # 去掉默认条
    # strip.text = element_text(face = "bold", size = 12), 
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  # 自定义背景填色：手动给分面上色
  # ggh4x::facet_wrap2(~ lineage, strip = ggh4x::strip_themed(
  #   background_x = ggh4x::elem_list_rect(fill = c("genome" = "white"))
  # )) + 
  coord_cartesian(ylim = c(-0.055, 0.105)) + 
  theme(panel.background = element_rect(fill = "#f9f9f9"))

cowplot::plot_grid(
  p1, p2,
  ncol = 2,               # 垂直排列（1 列）
  align = "h",            # 垂直对齐
  rel_widths = c(0.6, 1)   # 可调整两个图的高度比例
)

# ggsave('./plot/Lineage Prevalence/phylo_time_rates_of_gene_gain_and_loss_type2.pdf', width = 7, height = 4.5)



# dat <- rbind(
#   dat_soft %>% mutate(pangenome = 'soft_gene'), 
#   dat_rare %>% mutate(pangenome = 'rare_gene')
# ) %>%
#   pivot_longer(
#     cols = c(acc_gain, acc_loss),
#     names_to = "event",
#     values_to = "acc"
#   ) %>%
#   mutate(event = ifelse(event == "acc_gain", "gain", "loss"))
# 
# m <- stats::glm(acc ~ 0 + depth + core + lineage:core + lineage:depth + lineage:event:core + lineage:depth + event:depth + lineage:event:depth, data = dat %>% dplyr::filter(pangenome == "soft_gene"), family = "quasipoisson")
# 
# summary(m)
# 
# marginaleffects::slopes(m, variables = "core", by = c("lineage", "event")) %>% 
#   as.data.frame() %>% 
#   ggplot(aes(x = estimate, y = lineage, color = event)) +
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, position = position_dodge(width = 0.5)) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#   labs(
#     x = "margin effect",
#     y = "Lineage"
#   ) + 
#   coord_flip() + 
#   theme_minimal()
