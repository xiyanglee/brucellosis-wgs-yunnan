# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("lmc297/bactaxR"))
library(tidyverse)
library(ape)
library(phangorn)
library(ggplot2)
library(bactaxR)
library(ggdendro)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')



## 合并 fastANI 结果（由于文件数量过多，无法直接通配符合并）
# 
# find ./results/ -name "*.out" | while read FILE; do cat "$FILE" >> ./fastani_all_results.out; done
# 
## 从合并结果中去除完整路径，保留样本名
# 
# awk -F'\t' '{
#   split($1, a, "/"); sample1 = a[length(a)-1];
#   split($2, b, "/"); sample2 = b[length(b)-1];
#   print sample1 "\t" sample2 "\t" $3 "\t" $4 "\t" $5;
# }' fastani_all_results.out > fastani_results.out
# 
## 补全对称矩阵
# 
# awk -F '\t' '{
#   print $0;                                # 保留原始行
#   print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5; # 添加对称行：交换 sample1 和 sample2
# }' fastani_results.out > fastani_results_symmetric.out
# 
## 删除外类群
# 
# grep -v 'out_group' fastani_results_symmetric.out > fastani_results_symmetric.out.tmp
# mv fastani_results_symmetric.out.tmp fastani_results_symmetric.out



ani <- read.ANI(file = "./data/fastANI_result/fastani_results_symmetric.out")
summary(ani)

head(sort(table(ani@query[ani@ANI < 99.81]), decreasing = TRUE), 20)
# GCA_947242715.1 GCA_947242235.1 GCA_947241995.1 GCA_947242515.1 GCA_039649315.1 GCA_947242755.1 GCA_009823695.1 GCA_000250775.2 GCA_000250835.2 
#             816             771             749             748             571              81              12              11               9 

## 删除与其他样本 ANI 异常低的样本
remove_samples <- names(head(sort(table(ani@query[ani@ANI < 99.7]), decreasing = TRUE), 16))
keep_idx <- !(ani@query %in% remove_samples | ani@reference %in% remove_samples)

ani@ANI       <- ani@ANI[keep_idx]
ani@query     <- ani@query[keep_idx]
ani@reference <- ani@reference[keep_idx]

rm(keep_idx)


## 提取自函数 ANI.histogram
# hclust.method: Agglomeration method to pass to hclust. Defaults to "average".
# one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
as.dend <- function(bactaxRObject, hclust.method = "average") {
  fastani <- data.frame(bactaxRObject@query, bactaxRObject@reference, 
                        bactaxRObject@ANI)
  colnames(fastani) <- c("query", "reference", "ANI")
  s <- dcast(fastani, formula <- query ~ reference, value.var = "ANI")
  rownames(s) <- s$query
  s <- as.matrix(s[, !(colnames(s) == "query")])
  j <- matrix(data = 100, nrow = nrow(s), ncol = ncol(s))
  d <- j - s
  d.sym <- 0.5 * (d + t(d))
  d.dist <- as.dist(d.sym)
  h <- hclust(d = d.dist, method = hclust.method)
  dend <- as.dendrogram(h)
  return(dend)
}

dend <- as.dend(bactaxRObject = ani, hclust.method = "ward.D2")
phy  <- as.phylo(as.hclust(dend))


# 递归获取给定节点的所有父节点
getRecursivelyParentNode <- function(phy, node) {
  parents <- integer(0)
  current <- node
  while (!is.null(current <- phytools::getParent(phy, current))) {
    parents <- c(parents, current)
  }
  return(parents)
}

# 输入：phy 为 phylo 对象，tip.list 为目标分群的 tip 名称向量
getMrcaAndSubtreeNodes <- function(phy, tip.list) {
  # 获取 Most Recent Common Ancestors
  mrca_node  <- ape::getMRCA(phy, tip.list)
  # 获取所有父祖节点
  parent_nodes <- getRecursivelyParentNode(phy, mrca_node)
  # 获取所有子孙节点
  desc_nodes <- phytools::getDescendants(phy, mrca_node)
  desc_nodes <- desc_nodes[desc_nodes > phy$Nnode]
  all_nodes  <- sort(c(parent_nodes, mrca_node, desc_nodes))
  return(all_nodes)
}

getSubtreeNodes <- function(phy, mrca_node) {
  # 获取所有父祖节点
  parent_nodes <- getRecursivelyParentNode(phy, mrca_node)
  # 获取所有子孙节点
  desc_nodes <- phytools::getDescendants(phy, mrca_node)
  desc_nodes <- desc_nodes[desc_nodes > phy$Nnode]
  all_nodes  <- sort(c(parent_nodes, mrca_node, desc_nodes))
  return(all_nodes)
}

# ANI 数据：query x reference 的矩阵（从 bactaxRObject 提取）
getANIMatrix <- function(bactaxRObject) {
  fastani <- data.frame(bactaxRObject@query, bactaxRObject@reference, 
                        bactaxRObject@ANI)
  colnames(fastani) <- c("query", "reference", "ANI")
  s <- reshape2::dcast(fastani, query ~ reference, value.var = "ANI")
  rownames(s) <- s$query
  s <- s[, -1]
  s <- as.matrix(s)
  s[lower.tri(s)] <- t(s)[lower.tri(s)]  # 对称补全
  diag(s) <- 100
  return(s)
}

# 对一个节点，提取其子树下所有 tips，并计算 SD
computeNodeSD <- function(node, phy, ani_mat) {
  # phangorn::Descendants 获取不包含中间节点的子节点
  tips <- sort(phy$tip.label[unlist(phangorn::Descendants(phy, node, type = "tips"))])
  # 如果该节点下不足两个 tip（即无法构成 ANI 两两组合），返回 NA
  if (length(tips) < 2) return(NA_real_)
  #  生成所有两个 tip 的组合（两列矩阵，第一行为组合中的第一个样本名，第二行为第二个）
  comb <- combn(tips, 2)
  # map2_dbl() 是 purrr 包的函数：对两个等长向量 .x 和 .y，分别执行 .x[i], .y[i] → 某函数，返回 double 向量
  # 对每一对样本名从 ANI 矩阵中提取 ANI 值：ani_mat[.x, .y]
  # 得到当前子树中所有样本对之间的 ANI 值向量
  ani_vals <- purrr::map2_dbl(comb[1,], comb[2,], ~ ani_mat[.x, .y])
  # 计算并返回 ANI 值的标准差（即子树中遗传相似性的离散程度）
  return(sd(ani_vals, na.rm = TRUE))
}

computeNodeSD <- function(node, phy, ani_mat, method = c("iqr", "mad", "trimmed", "sd"), trim_prop = 0.1) {
  method <- match.arg(method)
  tips <- sort(phy$tip.label[unlist(phangorn::Descendants(phy, node, type = "tips"))])
  if (length(tips) < 2) return(NA_real_)
  
  comb <- combn(tips, 2)
  ani_vals <- purrr::map2_dbl(comb[1,], comb[2,], ~ ani_mat[.x, .y])
  
  # 去除 NA
  ani_vals <- ani_vals[!is.na(ani_vals)]
  if (length(ani_vals) < 2) return(NA_real_)
  
  # 使用不同方法计算稳健的“标准差”
  if (method == "iqr") {
    return(IQR(ani_vals) / 1.349)
  } else if (method == "mad") {
    return(mad(ani_vals, constant = 1.4826))  # 常数用于匹配正态分布 SD
  } else if (method == "trimmed") {
    lower <- quantile(ani_vals, probs = trim_prop)
    upper <- quantile(ani_vals, probs = 1 - trim_prop)
    trimmed_vals <- ani_vals[ani_vals >= lower & ani_vals <= upper]
    return(sd(trimmed_vals, na.rm = TRUE))
  } else {
    return(sd(ani_vals, na.rm = TRUE))  # 原始标准差
  }
}



# 使用 ggtree 可视化，并显示所有节点编号
# ggtree(phy) +
#   geom_text2(aes(label = node), hjust = -0.2, size = 3)


ani_mat <- getANIMatrix(ani)

## (1) 通过预先设置 tip list 获取 MRCA 节点
# tip.list <- partition.tree$tip[partition.tree$group.tree %in% c(4)]
# node_set <- getMrcaAndSubtreeNodes(phy, tip.list)

## (2) 直接定义 MRCA 节点
mrca_node <- 852
node_set <- getSubtreeNodes(phy, mrca_node)


sd_table <- tibble(
  node = node_set,
  sd = map_dbl(node_set, ~ computeNodeSD(.x, phy, ani_mat))
  # sd = map_dbl(node_set, ~ computeNodeSD(.x, phy, ani_mat, method = "mad"))
) %>% 
  mutate(node_type = case_when(
    node < {{ mrca_node }} ~ 'Parents', 
    node == {{ mrca_node }} ~ 'Current', 
    node > {{ mrca_node }} ~ 'Descendants'
  ))

if (nrow(sd_table) > 20) {
  sd_table <- bind_rows(
    sd_table[1:30, ], 
    sd_table[21:nrow(sd_table), ] %>% 
      summarise(
        node = 2000,
        sd = mean(sd, na.rm = TRUE),
        node_type = "Descendants"
      )
  )
}

sd_table %>% 
  mutate(node_name = paste0("n_", as.character(node + 16))) %>% 
  ggplot(aes(x = reorder(node_name, node), y = sd)) +
  geom_col(aes(fill = node_type)) + 
  geom_hline(yintercept = 0.02, linetype = "dashed", color = "black") + 
  labs(x = "Node index (DFS path)", y = "Standard deviation of ANI") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) + 
  coord_fixed(ratio = 100) + 
  scale_fill_manual(values = c("Parents" = "#3c4654", "Current" = "#108b96", "Descendants" = "#bdcbd2"), na.value = NA) + 
  theme(legend.position = "none")



################################################################################
## 循环绘图

# 定义 MRCA 节点和对应名称
mrca_named_list <- tibble(
  mrca_node = c(846, 840, 852, 845, 842, 844, 850, 851, 848),
  lineage_name = c('1.1', '1.2', '1.3', '2.1', '2.2', '2.3', '3.1', '3.2', '3.3')
)

# 创建输出文件夹（如果不存在）
dir.create("./plot/lineage_sd_plots", showWarnings = FALSE)

# 循环分析并保存图
for (i in seq_len(nrow(mrca_named_list))) {
  
  mrca_node <- mrca_named_list$mrca_node[i]
  lineage_name <- mrca_named_list$lineage_name[i]
  
  node_set <- getSubtreeNodes(phy, mrca_node)
  
  sd_table <- tibble(
    node = node_set,
    sd = map_dbl(node_set, ~ computeNodeSD(.x, phy, ani_mat))
  ) %>% 
    mutate(node_type = case_when(
      node < mrca_node ~ 'Parents', 
      node == mrca_node ~ 'Current', 
      node > mrca_node ~ 'Descendants'
    ))
  
  ## 30 for Lineage_xx.pdf
  # if (nrow(sd_table) > 30) {
  #   sd_table <- bind_rows(
  #     sd_table[1:30, ], 
  #     sd_table[31:nrow(sd_table), ] %>% 
  #       summarise(
  #         node = 3000,
  #         sd = mean(sd, na.rm = TRUE),
  #         node_type = "Descendants"
  #       )
  #   )
  # }
  ## 60 for Lineage_xx_deep.pdf
  if (nrow(sd_table) > 60) {
    sd_table <- bind_rows(
      sd_table[1:60, ],
      sd_table[61:nrow(sd_table), ] %>%
        summarise(
          node = 6000,
          sd = mean(sd, na.rm = TRUE),
          node_type = "Descendants"
        )
    )
  }
  
  # 绘图
  p <- sd_table %>% 
    mutate(node_name = paste0("n_", as.character(node + 16))) %>% 
    ggplot(aes(x = reorder(node_name, node), y = sd)) +
    geom_col(aes(fill = node_type)) + 
    geom_hline(yintercept = 0.02, linetype = "dashed", color = "black") + 
    labs(title = paste0("Lineage ", lineage_name),
         x = "Node index (DFS path)", 
         y = "Standard deviation of ANI") +
    theme_bw(base_size = 12) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) + 
    coord_fixed(ratio = 100) + 
    scale_fill_manual(values = c(
      "Parents" = "#3c4654", 
      "Current" = "#108b96", 
      "Descendants" = "#bdcbd2"
    ), na.value = NA) + 
    theme(legend.position = "none")
  
  # 保存图像
  ggsave(
    # filename = paste0("./plot/lineage_sd_plots/Lineage_", lineage_name, ".pdf"),
    # plot = p, height = 3, width = 8
    filename = paste0("./plot/lineage_sd_plots/Lineage_", lineage_name, "_deep.pdf"),
    plot = p, height = 3, width = 16
  )
}


# 使用 ggtree 可视化，并显示所有节点编号
ggtree(phy) + 
  geom_text2(aes(label = node, subset = node < 852), hjust = -0.2, size = 3)


ggtree(phy) + 
  # 添加内部节点编号，隐藏 tip 上的编号
  geom_text2(aes(label = paste0("n_", as.character(node + 16)), subset = node < 852 & node >= 835), 
             hjust = -0.2, size = 4, color = "#028760") + 
  theme_tree2() + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dotted"),
    panel.background = element_blank()
  ) + 
  # 添加 x 轴，显示距离差异度
  scale_x_continuous(name = "Pairwise ANI dissimilarity (%)", labels = ~ round(1 - ., 3))

# ggsave(filename = paste0("./plot/lineage_sd_plots/ANI_tree.pdf"), height = 8, width = 8)



################################################################################
## 查看具体节点和 BAPS 分类群关系（辅助）

load('./r_image/clustering_analysis_fastBAPS.RData')

# node ~ baps tree cluster ~ 'lineage'
# 846 ~ c(11, 12, 13, 14, 15) ~ '1.1'
# 840 ~ c(5, 6, 7, 8, 9) ~ '1.2'
# 852 ~ c(16, 10) ~ '1.3'
# 844 ~ 863 ~ '1.3a'
# 844 ~ 862 ~ '1.3b'
# 
# 845 ~ 17 ~ '2.1'
# 842 ~ 3 ~ '2.2'
# 844 ~ 18 ~ '2.3'
# 
# 850 ~ 1 ~ '3.1'
# 851 ~ 4 ~ '3.2'
# 848 ~ 2 ~ '3.3'

ggtree(raxml_tree)

# 使用 ggtree 可视化，并显示所有节点编号
ggtree(phy) +
  geom_text2(aes(label = node), hjust = -0.2, size = 3)

# 转换为 ggdendro 对象
dend_data <- dendro_data(dend)

dend_data$labels <- dend_data$labels %>% 
  mutate(
    highlight = case_match(
      label, 
      phy$tip.label[unlist(phangorn::Descendants(phy, 863, type = "tips"))] ~ 'labeld',
      # partition.tree$tip[partition.tree$group.tree %in% c(16, 10)] ~ 'labeld',
      .default = NA
    )
  )

# 绘图
ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey40") +
  geom_point(data = dend_data$labels,
             aes(x = x, y = y, color = highlight),
             size = 2) +
  coord_flip() +
  scale_y_reverse() +
  theme_minimal() +
  theme(legend.position = "none")  # 去掉图例
