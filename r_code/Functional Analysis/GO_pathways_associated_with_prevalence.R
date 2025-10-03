library(dplyr)
library(lme4)
library(car)
library(purrr)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')
# load('./r_image/enricher_db_panaroo.RData')
load('./r_image/enricher_db_bacteria_v2.RData')


prevalence_rate_by_year <- metadata.sample %>% 
  group_by(collection_year, lineage_level_2) %>% 
  summarise(covered_counties = n_distinct(county_gb), .groups = "drop")  %>%
  mutate(prevalence_rate = covered_counties / length(unique(metadata.sample$county_gb)))

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

# datmat.onehot <- read.table('./data/panaroo_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat.onehot <- read.table('./data/roary_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat.onehot <- normalize_datamat(datmat.onehot)
datmat.onehot <- datmat.onehot[grepl("^20", rownames(datmat.onehot)), ]


all_GO_terms <- unique(db.GO.emapper$TERM2GENE$TERM)

results <- map_dfr(all_GO_terms, function(go_term) {
  genes <- db.GO.emapper$TERM2GENE$GENE[db.GO.emapper$TERM2GENE$TERM == go_term]
  genes <- intersect(genes, colnames(datmat.onehot))
  if (length(genes) == 0) return(NULL)
  
  go_rate <- as.data.frame(rowMeans(datmat.onehot[, genes, drop = FALSE]))
  colnames(go_rate) <- 'go_rate'
  go_rate$assembly_accession <- rownames(go_rate)
  
  df <- left_join(metadata.sample, go_rate, by = 'assembly_accession') %>% 
    group_by(lineage_level_2, collection_year) %>%
    summarise(go_rate = mean(go_rate, na.rm = TRUE) * 100, .groups = 'drop') %>%
    left_join(., prevalence_rate_by_year, by = c('lineage_level_2', 'collection_year')) %>% 
    dplyr::rename(group = lineage_level_2) %>% 
    filter(!is.na(prevalence_rate), !is.na(go_rate))
  
  # 如果 go_rate 没有变异或数据不足，不建模
  if (nrow(df) < 5 || length(unique(df$group)) < 2 || var(df$go_rate) == 0) return(NULL)
  
  model <- tryCatch(lmer(prevalence_rate ~ go_rate + (1|group), data = df), error = function(e) return(NULL))
  if (is.null(model)) return(NULL)
  
  anova_res <- tryCatch(car::Anova(model), error = function(e) return(NULL))
  coef_res <- tryCatch(summary(model)$coefficients, error = function(e) return(NULL))
  
  if (is.null(anova_res) || is.null(coef_res)) return(NULL)
  if (!("go_rate" %in% rownames(anova_res)) || !("go_rate" %in% rownames(coef_res))) return(NULL)
  
  tibble(
    GO_term = go_term,
    Estimate = coef_res["go_rate", "Estimate"],
    Std_Error = coef_res["go_rate", "Std. Error"],
    p_value = anova_res["go_rate", "Pr(>Chisq)"]
  )
})

results$q_value <- p.adjust(results$p_value, method = "bonferroni")
results <- results %>% arrange(q_value)




plot_go_rate_vs_prevalence <- function(go_term, db.GO.emapper, datmat, metadata_yunnan, prevalence_rate_by_year) {
  # 提取基因列表
  genes <- db.GO.emapper$TERM2GENE$GENE[db.GO.emapper$TERM2GENE$TERM == go_term]
  genes <- intersect(genes, colnames(datmat))
  if (length(genes) == 0) {
    message("No genes found for ", go_term)
    return(NULL)
  }
  
  # 计算 go_rate
  go_rate <- as.data.frame(rowMeans(datmat[, genes, drop = FALSE]))
  colnames(go_rate) <- 'go_rate'
  go_rate$assembly_accession <- rownames(go_rate)
  
  # 合并数据
  df <- left_join(metadata.sample, go_rate, by = 'assembly_accession') %>% 
    group_by(lineage_level_2, collection_year) %>%
    summarise(go_rate = mean(go_rate, na.rm = TRUE), .groups = 'drop') %>%
    left_join(., prevalence_rate_by_year, by = c('lineage_level_2', 'collection_year')) %>% 
    dplyr::rename(group = lineage_level_2) %>% 
    filter(!is.na(prevalence_rate), !is.na(go_rate))
  
  if (nrow(df) < 5) {
    message("Not enough data points for ", go_term)
    return(NULL)
  }
  
  # 绘图
  p <- ggplot(df, aes(x = go_rate, y = prevalence_rate)) +
    geom_point(aes(color = group), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    labs(
      title = paste("Prevalence Rate vs GO Rate:", go_term),
      x = "Mean GO Pathway Abundance (go_rate)",
      y = "Prevalence Rate"
    ) +
    theme_minimal()
  
  return(p)
}

# plot_go_rate_vs_prevalence("GO:0006575", db.GO.emapper, datmat.onehot, metadata.sample, prevalence_rate_by_year)
# db.GO.emapper$TERM2NAME$NAME[db.GO.emapper$TERM2NAME$TERM == "GO:0006575"]
# db.GO.emapper$TERM2NAME$NAME[db.GO.emapper$TERM2NAME$TERM %in% results$GO_term[results$Estimate < 0 & results$q_value < 0.05]]


################################################################################
## 网络图 (风格 1)

plot_GO_network_WGCNA <- function(sig_terms, db.GO.emapper, results, datmat, edge_cutoff = 0.1) {
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(scales)
  
  # GO term → gene 映射
  term2gene <- db.GO.emapper$TERM2GENE %>%
    filter(TERM %in% sig_terms)
  
  valid_genes <- colnames(datmat)
  
  go_gene_list <- split(term2gene$GENE, term2gene$TERM) %>%
    map(~ intersect(.x, valid_genes)) %>%
    keep(~ length(.x) > 0)
  
  go_terms <- names(go_gene_list)
  
  # Jaccard 计算
  combn_df <- expand.grid(T1 = go_terms, T2 = go_terms, stringsAsFactors = FALSE) %>%
    filter(T1 < T2) %>%
    rowwise() %>%
    mutate(
      inter = length(intersect(go_gene_list[[T1]], go_gene_list[[T2]])),
      union = length(union(go_gene_list[[T1]], go_gene_list[[T2]])),
      jaccard = ifelse(union > 0, inter / union, 0)
    ) %>%
    ungroup() %>%
    filter(jaccard >= edge_cutoff)
  
  if (nrow(combn_df) == 0) {
    message("No edges with Jaccard ≥ ", edge_cutoff)
    return(NULL)
  }
  
  # 构建图
  g <- graph_from_data_frame(combn_df, directed = FALSE)
  
  # 节点属性：Estimate + 连接度
  V(g)$GO_term <- V(g)$name
  V(g)$Estimate <- results$Estimate[match(V(g)$GO_term, results$GO_term)]
  # V(g)$Estimate_plot <- ifelse(V(g)$Estimate > 50, 50, V(g)$Estimate)
  V(g)$Degree <- abs(V(g)$Estimate)  # degree(g), 可替换为 abs(Estimate) 作为 size 
  
  # WGCNA 风格图
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = jaccard, alpha = jaccard), color = "#383c3c") +
    geom_node_point(aes(size = Degree, color = Estimate)) +
    geom_node_text(aes(label = GO_term), repel = TRUE, size = 3, family = "sans") +
    scale_color_gradient2(low = muted("#165e83"), mid = "#fbfaf5", high = muted("#e83929"), midpoint = 0) +
    scale_edge_width(range = c(0.2, 2), guide = "none") +
    scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
    scale_size_continuous(range = c(3, 8)) +
    theme_void() +
    labs(
      # title = "GO Network (WGCNA Style)",
      color = "Effect Size",
      size = "Degree"
    )
}

plot_GO_network_WGCNA(
  sig_terms = results$GO_term[results$q_value < 0.05],
  db.GO.emapper = db.GO.emapper,
  results = results,
  datmat = datmat.onehot,
  edge_cutoff = 0.1
)


################################################################################
## 网络图 (风格 2)

plot_GO_network_WGCNA <- function(sig_terms, db.GO.emapper, results, datmat, edge_cutoff = 0.1, top_n_modules = 18) {
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(scales)
  library(tidygraph)
  
  # 自定义颜色
  col_g <- "#C1C1C1"
  cols <- c("#DEB99B", "#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D",
            "#6E4821", "#A4B423", "#C094DF", "#DC95D8", "#326530", "#50C0C9",
            "#67C021", "#DC69AF", "#8C384F", "#30455C", "#F96C72", "#5ED2BF")
  
  # GO term → gene 映射
  term2gene <- db.GO.emapper$TERM2GENE %>%
    filter(TERM %in% sig_terms)
  
  valid_genes <- colnames(datmat)
  go_gene_list <- split(term2gene$GENE, term2gene$TERM) %>%
    map(~ intersect(.x, valid_genes)) %>%
    keep(~ length(.x) > 0)
  go_terms <- names(go_gene_list)
  
  # Jaccard 相似度计算
  combn_df <- expand.grid(T1 = go_terms, T2 = go_terms, stringsAsFactors = FALSE) %>%
    filter(T1 < T2) %>%
    rowwise() %>%
    mutate(
      inter = length(intersect(go_gene_list[[T1]], go_gene_list[[T2]])),
      union = length(union(go_gene_list[[T1]], go_gene_list[[T2]])),
      jaccard = ifelse(union > 0, inter / union, 0)
    ) %>%
    ungroup() %>%
    filter(jaccard >= edge_cutoff)
  
  if (nrow(combn_df) == 0) {
    message("No edges with Jaccard ≥ ", edge_cutoff)
    return(NULL)
  }
  
  # 构建图
  g <- graph_from_data_frame(combn_df, directed = FALSE)
  V(g)$GO_term <- V(g)$name
  V(g)$Estimate <- results$Estimate[match(V(g)$GO_term, results$GO_term)]
  V(g)$Degree <- abs(V(g)$Estimate)
  
  # 模块识别
  set.seed(7)
  V(g)$modularity <- membership(cluster_fast_greedy(g))
  
  # 模块颜色设定
  modu_sort <- sort(table(V(g)$modularity), decreasing = TRUE)
  modu_name <- names(modu_sort)[1:min(top_n_modules, length(modu_sort))]
  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  
  V(g)$color <- ifelse(
    V(g)$modularity %in% modu_name,
    modu_cols[as.character(V(g)$modularity)],
    col_g
  )
  
  # 边颜色设定
  edge_df <- as_data_frame(g, what = "edges")
  edge_df$mod1 <- V(g)$modularity[match(edge_df$from, V(g)$name)]
  edge_df$mod2 <- V(g)$modularity[match(edge_df$to, V(g)$name)]
  edge_df$color <- ifelse(
    edge_df$mod1 == edge_df$mod2 & edge_df$mod1 %in% modu_name,
    modu_cols[as.character(edge_df$mod1)],
    col_g
  )
  E(g)$color <- edge_df$color
  E(g)$jaccard <- edge_df$jaccard
  
  # 计算布局
  layout_coords <- layout_with_fr(g, niter = 999, grid = "nogrid")
  
  # 转换为 tbl_graph 并添加节点属性（Degree, GO_term, color, modularity）
  g_tbl <- as_tbl_graph(g) %>%
    mutate(
      x = layout_coords[, 1],
      y = layout_coords[, 2],
      Degree = V(g)$Degree,
      GO_term = V(g)$GO_term,
      color = V(g)$color,
      modularity = V(g)$modularity
    )
  
  # 使用 create_layout 显式生成 layout
  g_layout <- create_layout(g_tbl, layout = "manual", x = x, y = y)
  
  # 绘图
  ggraph(g_layout) +
    geom_edge_link(aes(edge_alpha = jaccard, edge_width = jaccard, edge_colour = I(color))) +
    geom_node_point(aes(color = I(color), size = Degree)) +
    geom_node_text(aes(label = GO_term), repel = TRUE, size = 3, family = "sans") +
    scale_edge_width(range = c(0.2, 2), guide = "none") +
    scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
    scale_size_continuous(range = c(3, 8)) +
    theme_void() +
    labs(
      color = "Module",
      size = "Effect Size"
    )
}

plot_GO_network_WGCNA(
  sig_terms = results$GO_term[results$q_value < 0.05],
  db.GO.emapper = db.GO.emapper,
  results = results,
  datmat = datmat.onehot,
  edge_cutoff = 0.1
)

################################################################################
## 网络图 (风格 3)

plot_GO_network_WGCNA <- function(sig_terms, db.GO.emapper, results, datmat, 
                                  edge_cutoff = 0.1, top_n_modules = 18,
                                  size.range = c(3, 8), edge.range = c(0.2, 2)) {
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(scales)
  library(tidygraph)
  library(ggalt)
  library(rlang)
  
  # 自定义颜色
  col_g <- "#203744"
  cols <- c("#0f2350", "#206864", "#a22041", "#55295b", "#164a84")
  
  # GO term → gene 映射
  term2gene <- db.GO.emapper$TERM2GENE %>%
    filter(TERM %in% sig_terms)
  valid_genes <- colnames(datmat)
  go_gene_list <- split(term2gene$GENE, term2gene$TERM) %>%
    map(~ intersect(.x, valid_genes)) %>%
    keep(~ length(.x) > 0)
  go_terms <- names(go_gene_list)
  
  # Jaccard 相似度计算
  combn_df <- expand.grid(T1 = go_terms, T2 = go_terms, stringsAsFactors = FALSE) %>%
    filter(T1 < T2) %>%
    rowwise() %>%
    mutate(
      inter = length(intersect(go_gene_list[[T1]], go_gene_list[[T2]])),
      union = length(union(go_gene_list[[T1]], go_gene_list[[T2]])),
      jaccard = ifelse(union > 0, inter / union, 0)
    ) %>%
    ungroup() %>%
    filter(jaccard >= edge_cutoff)
  
  if (nrow(combn_df) == 0) {
    message("No edges with Jaccard ≥ ", edge_cutoff)
    return(NULL)
  }
  
  # 构建图
  g <- graph_from_data_frame(combn_df, directed = FALSE)
  V(g)$GO_term <- V(g)$name
  V(g)$Estimate <- results$Estimate[match(V(g)$GO_term, results$GO_term)]
  V(g)$Degree <- abs(V(g)$Estimate)
  
  # 模块识别
  set.seed(7)
  V(g)$modularity <- membership(cluster_fast_greedy(g))
  modu_sort <- sort(table(V(g)$modularity), decreasing = TRUE)
  modu_name <- names(modu_sort)[1:min(top_n_modules, length(modu_sort))]
  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  
  # 节点颜色与描边
  V(g)$modu_color <- ifelse(
    V(g)$modularity %in% modu_name,
    modu_cols[as.character(V(g)$modularity)],
    col_g
  )
  
  # 边属性：同模块上色，否则灰色
  edge_df <- as_data_frame(g, what = "edges")
  mod1 <- V(g)$modularity[match(edge_df$from, V(g)$name)]
  mod2 <- V(g)$modularity[match(edge_df$to, V(g)$name)]
  same_mod <- mod1 == mod2 & mod1 %in% as.integer(modu_name)
  edge_df$color <- ifelse(same_mod, modu_cols[as.character(mod1)], col_g)
  E(g)$jaccard <- edge_df$jaccard
  E(g)$color <- edge_df$color
  
  # 计算布局
  # layout_coords <- layout_with_fr(g, niter = 999, grid = "nogrid")
  # layout_coords <- layout_with_fr(g, niter = 9999)
  layout_coords <- layout_with_fr(g, weights = E(g)$jaccard ** 1.8, niter = 999)
  
  # 构建 tidygraph
  g_tbl <- as_tbl_graph(g) %>%
    mutate(
      x = layout_coords[, 1],
      y = layout_coords[, 2],
      Degree = V(g)$Degree,
      GO_term = V(g)$GO_term,
      Estimate = V(g)$Estimate,
      Estimate_plot = pmax(pmin(V(g)$Estimate, 2), -2),  # 将 >2 设为2，<−2 设为−2
      modularity = V(g)$modularity,
      modu_color = V(g)$modu_color
    )
  
  # 创建 layout
  g_layout <- create_layout(g_tbl, layout = "manual", x = x, y = y)
  
  # 获取模块坐标
  mod_groups <- as.data.frame(g_layout) %>%
    filter(modularity %in% as.numeric(modu_name))
  
  # 主图绘制
  p <- ggraph(g_layout) +
    # 模块边界
    purrr::map(unique(mod_groups$modularity), function(mod_id) {
      mod_df <- mod_groups[mod_groups$modularity == mod_id, ]
      ggalt::geom_encircle(
        data = mod_df,
        aes(x = x, y = y),
        color = modu_cols[as.character(mod_id)],
        # s_shape = 1,
        expand = 0.05,
        alpha = 1,
        size = 2, 
        linetype = "dashed",
        inherit.aes = FALSE
      )
    }) +
    # geom_edge_link(aes(edge_width = jaccard, edge_colour = I(color)), alpha = 1) +
    geom_edge_link(aes(edge_width = jaccard), edge_colour = "black", alpha = 1) +
    # geom_node_point(aes(size = Degree, fill = Estimate, color = I(modu_color)), shape = 21, stroke = 1) + 
    geom_node_point(aes(size = Degree, fill = Estimate_plot), color = "black", shape = 21, stroke = 1) +
    # geom_node_text(aes(label = GO_term), repel = TRUE, size = 3, family = "sans") +
    # scale_fill_gradient2(low = muted("#4c6473"), mid = "#fbfaf5", high = muted("#d7003a"), midpoint = 0) + 
    scale_fill_gradient2(low = "#12507b", mid = "#efefef", high = "#c12c1f", midpoint = 0) + 
    scale_edge_width(range = edge.range, guide = "none") +
    scale_size_continuous(range = size.range) +
    theme_void() +
    labs(
      fill = "Estimate",
      size = "Effect Size",
      color = "Module Border"
    )
  
  return(p)
}

plot_GO_network_WGCNA(
  sig_terms = results$GO_term[results$q_value < 0.05],
  db.GO.emapper = db.GO.emapper,
  results = results,
  datmat = datmat.onehot,
  edge_cutoff = 0.1, 
  size.range = c(2, 5), 
  edge.range = c(0.1, 1)
)

################################################################################
## 网络图 (风格 3, 并返回 module, 最终版)

plot_GO_network_WGCNA <- function(sig_terms, db.GO.emapper, results, datmat, 
                                  edge_cutoff = 0.1, top_n_modules = 18,
                                  size.range = c(3, 8), edge.range = c(0.2, 2)) {
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(scales)
  library(tidygraph)
  library(ggalt)
  library(rlang)
  
  # 自定义颜色
  col_g <- "#203744"
  cols <- c("#0f2350", "#206864", "#a22041", "#55295b", "#164a84")
  
  # GO term → gene 映射
  term2gene <- db.GO.emapper$TERM2GENE %>%
    filter(TERM %in% sig_terms)
  valid_genes <- colnames(datmat)
  go_gene_list <- split(term2gene$GENE, term2gene$TERM) %>%
    map(~ intersect(.x, valid_genes)) %>%
    keep(~ length(.x) > 0)
  go_terms <- names(go_gene_list)
  
  # Jaccard 相似度计算
  combn_df <- expand.grid(T1 = go_terms, T2 = go_terms, stringsAsFactors = FALSE) %>%
    filter(T1 < T2) %>%
    rowwise() %>%
    mutate(
      inter = length(intersect(go_gene_list[[T1]], go_gene_list[[T2]])),
      union = length(union(go_gene_list[[T1]], go_gene_list[[T2]])),
      jaccard = ifelse(union > 0, inter / union, 0)
    ) %>%
    ungroup() %>%
    filter(jaccard >= edge_cutoff)
  
  if (nrow(combn_df) == 0) {
    message("No edges with Jaccard ≥ ", edge_cutoff)
    return(NULL)
  }
  
  # 构建图
  g <- graph_from_data_frame(combn_df, directed = FALSE)
  V(g)$GO_term <- V(g)$name
  V(g)$Estimate <- results$Estimate[match(V(g)$GO_term, results$GO_term)]
  V(g)$Degree <- abs(V(g)$Estimate)
  
  # 模块识别
  set.seed(7)
  V(g)$modularity <- membership(cluster_fast_greedy(g))
  modu_sort <- sort(table(V(g)$modularity), decreasing = TRUE)
  modu_name <- names(modu_sort)[1:min(top_n_modules, length(modu_sort))]
  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  
  # 节点颜色与描边
  V(g)$modu_color <- ifelse(
    V(g)$modularity %in% modu_name,
    modu_cols[as.character(V(g)$modularity)],
    col_g
  )
  
  # 边属性
  edge_df <- as_data_frame(g, what = "edges")
  mod1 <- V(g)$modularity[match(edge_df$from, V(g)$name)]
  mod2 <- V(g)$modularity[match(edge_df$to, V(g)$name)]
  same_mod <- mod1 == mod2 & mod1 %in% as.integer(modu_name)
  edge_df$color <- ifelse(same_mod, modu_cols[as.character(mod1)], col_g)
  E(g)$jaccard <- edge_df$jaccard
  E(g)$color <- edge_df$color
  
  # 计算布局（用 jaccard 控制布局）
  layout_coords <- layout_with_fr(g, weights = E(g)$jaccard ^ 1.8, niter = 999)
  
  # 构建 tidygraph 对象
  g_tbl <- as_tbl_graph(g) %>%
    mutate(
      x = layout_coords[, 1],
      y = layout_coords[, 2],
      Degree = V(g)$Degree,
      GO_term = V(g)$GO_term,
      Estimate = V(g)$Estimate,
      Estimate_plot = pmax(pmin(V(g)$Estimate, 2), -2),
      modularity = V(g)$modularity,
      modu_color = V(g)$modu_color
    )
  
  # 创建 layout
  g_layout <- create_layout(g_tbl, layout = "manual", x = x, y = y)
  
  # 模块坐标与标签
  mod_groups <- as.data.frame(g_layout) %>%
    filter(modularity %in% as.numeric(modu_name))
  
  label_df <- mod_groups %>%
    group_by(modularity) %>%
    summarise(x = median(x), y = median(y)) %>%
    mutate(label = paste0("Module ", modularity))
  
  # 绘图
  p <- ggraph(g_layout) +
    purrr::map(unique(mod_groups$modularity), function(mod_id) {
      mod_df <- mod_groups[mod_groups$modularity == mod_id, ]
      ggalt::geom_encircle(
        data = mod_df,
        aes(x = x, y = y),
        color = modu_cols[as.character(mod_id)],
        expand = 0.05,
        alpha = 1,
        size = 2, 
        linetype = "dashed",
        inherit.aes = FALSE
      )
    }) +
    geom_edge_link(aes(edge_width = jaccard, edge_color = jaccard, edge_alpha = jaccard)) + # , alpha = 1
    geom_node_point(aes(size = Degree, fill = Estimate_plot), color = "#2c2f3b", shape = 21, stroke = 1) +
    geom_text(data = label_df, aes(x = x, y = y, label = label), size = 4, fontface = "bold", color = "black") + 
    scale_edge_color_continuous(
      name = "Jaccard Index",
      low = "grey95", high = "#151d29",
      limits = c(0, 1)
    ) + 
    scale_fill_gradient2(
      name = "Effect size", 
      low = "#2792c3", mid = "#efefef", high = "#e9546b", 
      midpoint = 0
    ) + 
    # scale_edge_width(range = edge.range, guide = "none") + 
    scale_edge_alpha(
      range = c(0.25, 1),
      guide = "none"
    ) +
    scale_edge_width(
      name = "Jaccard Index",
      range = edge.range,
      breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)
    ) +
    scale_size_continuous(
      name = "Effect Size",
      range = size.range,
      breaks = c(0.1, 0.5, 1, 2, 3)
    ) + 
    theme_void() +
    labs(
      fill = "Estimate",
      # size = "Effect Size",
      color = "Module Border"
    )
  
  # 返回模块-GO term 表
  module_go_terms <- as.data.frame(g_layout) %>%
    filter(modularity %in% as.numeric(modu_name)) %>%
    group_by(modularity) %>%
    summarise(go_terms = list(GO_term)) %>%
    arrange(modularity)
  
  return(list(
    plot = p,
    module_go_terms = module_go_terms
  ))
}

res <- plot_GO_network_WGCNA(
  sig_terms = results$GO_term[results$q_value < 0.05],
  db.GO.emapper = db.GO.emapper,
  results = results,
  datmat = datmat.onehot,
  edge_cutoff = 0.1, 
  size.range = c(2, 5), 
  edge.range = c(0.5, 1)
)

p <- res$plot                       # 查看绘图
module_info <- res$module_go_terms  # 查看每个模块对应 GO terms
# rm(res)


module_info$go_terms_named <- purrr::map_chr(
  module_info$go_terms,
  function(term_list) {
    # 过滤出当前模块的TERM及其NAME
    term_names <- db.GO.emapper$TERM2NAME %>%
      filter(TERM %in% term_list) %>%
      mutate(label = paste0(NAME, " (", TERM, ")")) %>%
      pull(label)
    
    # 合并为一行字符串
    paste(term_names, collapse = "; ")
  }
)

module_info$go_terms <- NULL

# ggsave('./plot/Lineage Prevalence/GO_network.pdf', p, width = 5.5, height = 4)
# write.csv(module_info, './result/Lineage Prevalence/GO_module_info.csv', row.names = FALSE)
# save.image('./r_image/GO_pathways_associated_with_prevalence.RData')
