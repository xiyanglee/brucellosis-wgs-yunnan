# withr::with_libpaths("/public/r_share_library/4.1", remotes::install_github("gtonkinhill/panstripe", build_vignettes = TRUE))
library(panstripe)
library(ape)

setwd('/Users/xiyangli/Lab/Project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')
# setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

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

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.2"], tree$tip.label)
dat_l12 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.2")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.3"], tree$tip.label)
dat_l13 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.3")


dat_genome <- rbind(dat_l11, dat_l12, dat_l13)
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

# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0030556  0.0003817   8.006 1.03e-13 ***
# core            0.0016966  0.0032674   0.519  0.60418    
# core:lineage1.1 0.0122996  0.0038523   3.193  0.00164 ** 
# core:lineage1.2 0.0535296  0.1119878   0.478  0.63319    
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for quasipoisson family taken to be 44.73879)

#     Null deviance: 5787.9  on 200  degrees of freedom
# Residual deviance: 3163.5  on 196  degrees of freedom
# AIC: NA

# Number of Fisher Scoring iterations: 7
summary(m_loss)
# Call:
# stats::glm(formula = acc_loss ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_genome)

# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0029513  0.0003971   7.432 3.24e-12 ***
# core            0.0017784  0.0031969   0.556  0.57864    
# core:lineage1.1 0.0121478  0.0038101   3.188  0.00167 ** 
# core:lineage1.2 0.0458322  0.1216681   0.377  0.70681    
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for quasipoisson family taken to be 44.96909)

#     Null deviance: 5506.5  on 200  degrees of freedom
# Residual deviance: 3107.4  on 196  degrees of freedom
# AIC: NA

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
#   term lineage   estimate  std.error statistic     p.value   s.value    conf.low  conf.high      events gene_type
# 1 core     1.1 0.05618125 0.01755132 3.2009693 0.001369661 9.5119649  0.02178129 0.09058121 gain events    genome
# 2 core     1.2 0.10662376 0.36839077 0.2894311 0.772251467 0.3728574 -0.61540887 0.82865639 gain events    genome
# 3 core     1.3 0.01991440 0.03979980 0.5003643 0.616818578 0.6970819 -0.05809177 0.09792057 gain events    genome
# 4 core     1.1 0.05425955 0.01729450 3.1373882 0.001704603 9.1963488  0.02036296 0.08815614 loss events    genome
# 5 core     1.2 0.08346226 0.33758489 0.2472334 0.804727608 0.3134276 -0.57819197 0.74511648 loss events    genome
# 6 core     1.3 0.01926031 0.03618674 0.5322476 0.594554516 0.7501190 -0.05166440 0.09018502 loss events    genome


################################################################################
## -> soft-gene
pa <- datmat.onehot[, softgene_names]
tree <- beast_tree@phylo

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.1"], tree$tip.label)
dat_l11 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.1")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.2"], tree$tip.label)
dat_l12 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.2")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.3"], tree$tip.label)
dat_l13 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.3")


dat_soft <- rbind(dat_l11, dat_l12, dat_l13)
dat_soft$lineage <- as.factor(dat_soft$lineage)

m_gain <- stats::glm(acc_gain ~ 0 + depth + core + lineage:core, data = dat_soft, family = "quasipoisson")
m_loss <- stats::glm(acc_loss ~ 0 + depth + core + lineage:core, data = dat_soft, family = "quasipoisson")
summary(m_gain)
# Call:
# stats::glm(formula = acc_gain ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_soft)

# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0027215  0.0003595   7.570 1.42e-12 ***
# core            0.0015446  0.0028875   0.535  0.59332    
# core:lineage1.1 0.0104059  0.0036556   2.847  0.00489 ** 
# core:lineage1.2 0.0372115  0.1104592   0.337  0.73657    
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for quasipoisson family taken to be 30.80328)

#     Null deviance: 3961.5  on 200  degrees of freedom
# Residual deviance: 2396.5  on 196  degrees of freedom
# AIC: NA

# Number of Fisher Scoring iterations: 7
summary(m_loss)
# Call:
# stats::glm(formula = acc_loss ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_soft)

# Coefficients: (1 not defined because of singularities)
#                  Estimate Std. Error t value Pr(>|t|)    
# depth           0.0026790  0.0003908   6.855 9.06e-11 ***
# core            0.0016163  0.0030303   0.533   0.5944    
# core:lineage1.1 0.0118125  0.0036210   3.262   0.0013 ** 
# core:lineage1.2 0.0317146  0.1254817   0.253   0.8007    
# core:lineage1.3        NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for quasipoisson family taken to be 35.54695)

#     Null deviance: 4440.5  on 200  degrees of freedom
# Residual deviance: 2630.8  on 196  degrees of freedom
# AIC: NA

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
#   term lineage   estimate  std.error statistic     p.value   s.value    conf.low  conf.high      events gene_type
# 1 core     1.1 0.03852926 0.01273984 3.0243125 0.002491989 8.6484867  0.01355963 0.06349889 gain events soft_gene
# 2 core     1.2 0.06092295 0.25222498 0.2415421 0.809135006 0.3055477 -0.43342892 0.55527482 gain events soft_gene
# 3 core     1.3 0.01383512 0.02685420 0.5151938 0.606417609 0.7216164 -0.03879814 0.06646837 gain events soft_gene
# 4 core     1.1 0.04727101 0.01444748 3.2719208 0.001068195 9.8706092  0.01895447 0.07558755 loss events soft_gene
# 5 core     1.2 0.04919149 0.25518115 0.1927709 0.847138424 0.2393304 -0.45095436 0.54933734 loss events soft_gene
# 6 core     1.3 0.01402870 0.02744911 0.5110802 0.609294869 0.7147875 -0.03977057 0.06782796 loss events soft_gene


################################################################################
## -> rare-gene
pa <- datmat.onehot[, raregene_names]
tree <- beast_tree@phylo

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.1"], tree$tip.label)
dat_l11 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.1")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.2"], tree$tip.label)
dat_l12 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.2")

tips_to_keep <- intersect(metadata.sample$assembly_accession[metadata.sample$lineage_level_2 == "1.3"], tree$tip.label)
dat_l13 <- get_sub_fit_data(pa, tree, tips_to_keep, use_root_depth = FALSE) %>% mutate(lineage = "1.3")

dat_rare <- rbind(dat_l11, dat_l12, dat_l13)
dat_rare$lineage <- as.factor(dat_rare$lineage)

m_gain <- stats::glm(acc_gain ~ 0 + depth + core + lineage:core, data = dat_rare, family = "quasipoisson")
m_loss <- stats::glm(acc_loss ~ 0 + depth + core + lineage:core, data = dat_rare, family = "quasipoisson")

summary(m_gain)
# Call:
# stats::glm(formula = acc_gain ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_rare)

# Coefficients: (1 not defined because of singularities)
#                   Estimate Std. Error t value Pr(>|t|)
# depth            8.721e-05  6.082e-04   0.143    0.886
# core            -2.264e-03  7.673e-03  -0.295    0.768
# core:lineage1.1  1.292e-02  7.881e-03   1.640    0.103
# core:lineage1.2 -2.894e-02  1.477e-01  -0.196    0.845
# core:lineage1.3         NA         NA      NA       NA

# (Dispersion parameter for quasipoisson family taken to be 14.04019)

#     Null deviance: 1206.62  on 200  degrees of freedom
# Residual deviance:  993.74  on 196  degrees of freedom
# AIC: NA

# Number of Fisher Scoring iterations: 7
summary(m_loss)
# Call:
# stats::glm(formula = acc_loss ~ 0 + depth + core + lineage:core, 
#     family = "quasipoisson", data = dat_rare)

# Coefficients: (1 not defined because of singularities)
#                   Estimate Std. Error t value Pr(>|t|)
# depth           -0.0001196  0.0006474  -0.185    0.854
# core            -0.0005957  0.0043128  -0.138    0.890
# core:lineage1.1  0.0039706  0.0062930   0.631    0.529
# core:lineage1.2 -0.0533276  0.1762463  -0.303    0.763
# core:lineage1.3         NA         NA      NA       NA

# (Dispersion parameter for quasipoisson family taken to be 13.54359)

#     Null deviance: 865.57  on 200  degrees of freedom
# Residual deviance: 857.58  on 196  degrees of freedom
# AIC: NA

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
#   term lineage      estimate   std.error  statistic     p.value   s.value     conf.low   conf.high      events gene_type
# 1 core     1.1  0.0166425366 0.005589551  2.9774373 0.002906691 8.4264066  0.005687219 0.027597855 gain events rare_gene
# 2 core     1.2 -0.0236656467 0.083057443 -0.2849311 0.775696945 0.3664350 -0.186455243 0.139123949 gain events rare_gene
# 3 core     1.3 -0.0022305693 0.007382931 -0.3021252 0.762556649 0.3910836 -0.016700849 0.012239710 gain events rare_gene
# 4 core     1.1  0.0035521930 0.005306065  0.6694590 0.503202747 0.9907883 -0.006847504 0.013951890 loss events rare_gene
# 5 core     1.2 -0.0340854754 0.066473596 -0.5127671 0.608114221 0.7175858 -0.164371329 0.096200378 loss events rare_gene
# 6 core     1.3 -0.0005266832 0.003737846 -0.1409056 0.887944551 0.1714585 -0.007852726 0.006799359 loss events rare_gene



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
  theme_classic(base_size = 14) + 
  coord_cartesian(ylim = c(-0.05, 0.12))

# ggsave('./plot/Lineage Prevalence/phylo_time_rates_of_gene_gain_and_loss.pdf', width = 9, height = 3.1)

# save.image('./r_image/phylo_time_rates_of_gene_gain_and_loss.RData')
# load('./r_image/phylo_time_rates_of_gene_gain_and_loss.RData')

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
    # strip.background = element_blank(),  # 去掉默认条
    strip.text = element_text(color = "white", face = "bold", size = 12), 
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  # 自定义背景填色：手动给分面上色
  # ggh4x::facet_wrap2(~ lineage, strip = ggh4x::strip_themed(
  #   background_x = ggh4x::elem_list_rect(fill = c("1.1" = "#7EBB8F", "1.3" = "#886A9E"))
  # )) +
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

# ggsave('./plot/Lineage Prevalence/phylo_time_rates_of_gene_gain_and_loss_type3.pdf', width = 9, height = 4.5)

tmp_plot <- rbind(
  soft_gene_margin_effect, 
  rare_gene_margin_effect
) %>% 
  mutate(
    gene_type = factor(gene_type, levels = c("genome", "soft_gene", "rare_gene")), 
    sig_group = ifelse(p.value < 0.05, "p < 0.05", "p ≥ 0.05"),
    sig_group = factor(sig_group, levels = c("p < 0.05", "p ≥ 0.05"))
  ) %>% 
  dplyr::filter(lineage != "1.2")

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
    # strip.background = element_blank(),  # 去掉默认条
    strip.text = element_text(color = "white", face = "bold", size = 12), 
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  # 自定义背景填色：手动给分面上色
  # ggh4x::facet_wrap2(~ lineage, strip = ggh4x::strip_themed(
  #   background_x = ggh4x::elem_list_rect(fill = c("1.1" = "#7EBB8F", "1.3" = "#886A9E"))
  # )) +
  coord_cartesian(ylim = c(-0.055, 0.105)) + 
  theme(panel.background = element_rect(fill = "#f9f9f9"))

tmp_plot <- genome_margin_effect %>% 
  mutate(
    gene_type = factor(gene_type, levels = c("genome", "soft_gene", "rare_gene")), 
    sig_group = ifelse(p.value < 0.05, "p < 0.05", "p ≥ 0.05"),
    sig_group = factor(sig_group, levels = c("p < 0.05", "p ≥ 0.05"))
  ) %>% 
  dplyr::filter(lineage != "1.2")

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

rbind(
  genome_margin_effect, 
  soft_gene_margin_effect, 
  rare_gene_margin_effect
) %>% 
  mutate(
    gene_type = factor(gene_type, levels = c("genome", "soft_gene", "rare_gene")), 
    sig_group = ifelse(p.value < 0.05, "p < 0.05", "p ≥ 0.05"),
    sig_group = factor(sig_group, levels = c("p < 0.05", "p ≥ 0.05"))
  ) %>% 
  write.csv(., "./r_plot_data/figure_3d&e~margin_effect.csv", row.names = FALSE)


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
