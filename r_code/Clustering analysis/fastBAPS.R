# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("lmc297/bactaxR"))
library(fastbaps)
library(ape)
library(Matrix)
library(tidyverse)
library(ggtree)
library(dendextend) # trans tree to as.dendrogram

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')


## 通过 gubbins 输出的进化树检验异常枝
## outgroup:      'GCA_000369945.1'
## 异常分类群:    'GCA_000250835.2', 'GCA_009823695.1', 'GCA_000250775.2'
## 异常枝 (严重): 'GCA_000369945.1'

# tree <- ape::read.tree('./data/coregene_merge_snippy_gubbins/gubbins.final_tree.tre')
# tree <- root(tree, outgroup = 'GCA_000369945.1', resolve.root = TRUE)
# tree <- drop.tip(tree, c('GCA_000369945.1'))
# tree <- drop.tip(tree, c('GCA_000250835.2', 'GCA_009823695.1', 'GCA_000250775.2'))
# ggtree(tree)


## 读取进化树
raxml_tree <- ape::read.tree('./data/raxml_gubbins_contig/RAxML_bestTree.raxml_tree')
raxml_tree <- root(raxml_tree, outgroup = 'GCA_000369945.1', resolve.root = TRUE)
raxml_tree <- drop.tip(raxml_tree, c('GCA_000369945.1'))

## 将 Reference 替换回 GCA_000007125.1
"Reference" %in% raxml_tree$tip.label
# [1] TRUE
raxml_tree$tip.label[raxml_tree$tip.label == "Reference"] <- "GCA_000007125.1"
"GCA_000007125.1" %in% raxml_tree$tip.label
# [1] TRUE

## 修剪异常枝
# 获取所有 tip 的分支长度
lengths.tip_branch <- sapply(raxml_tree$tip.label, function(N) {
  target_index <- which(raxml_tree$tip.label == N)
  edge_index <- which(raxml_tree$edge[, 2] == target_index)
  raxml_tree$edge.length[edge_index]
})
head(sort(lengths.tip_branch, decreasing = TRUE), 10)
# GCA_039649315.1 GCA_039649395.1 GCA_900446075.1 GCA_039747915.1 GCA_900236575.1 GCA_000370745.1 GCA_000370785.1 GCA_039649475.1 GCA_002191415.1 GCA_000370825.1 
#     0.133685006     0.017536783     0.014191229     0.011050348     0.010545138     0.009347372     0.008681874     0.006586379     0.006471963     0.005482123 

tips.length_edit <- names(sort(lengths.tip_branch, decreasing = TRUE))[1:5]
tips.length_edit <- setdiff(tips.length_edit, 'GCA_900446075.1')
new_length <- unname(sort(lengths.tip_branch, decreasing = TRUE)[6])
for (N in tips.length_edit) {
  raxml_tree$edge.length[which(raxml_tree$edge[, 2] == which(raxml_tree$tip.label == N))] <- new_length
}; rm(N, new_length)

ggtree(raxml_tree)
ggtree(ape::compute.brlen(raxml_tree, method = "Grafen"))
ggtree(raxml_tree, layout = "ape")



sparse.base <- import_fasta_sparse_nt('./data/snippy_contig_gubbins/clean.core.aln')
sparse.base[["snp.matrix"]]@Dimnames[[2]][sparse.base[["snp.matrix"]]@Dimnames[[2]] == "Reference"] <- "GCA_000007125.1"

## 首先进行初步聚类，之后使用 best_baps_partition 组合较小的 clusters， 以最大化 BAPS likelihood
##  初步聚类可以通过以下方法

## (1) 基于 SNP  距离 (ref: ?best_baps_partition)
# d <- snp_dist(sparse.data)
# d <- as.dist(d/max(d))
# h <- hclust(d, method="ward.D2")

## (2) 通过进化树信息 (ref: ?best_baps_partition)
# iqtree <- phytools::read.newick(...)
# h <- phytools::midpoint.root(iqtree)

## (3) 使用 fast_baps 聚类 (ref: https://github.com/dvorikus/Microcoleus-population-genomics/blob/main/Clustering%20analysis/fastBAPS.R)
# baps.hc <- fast_baps(sparse.data, k.init=80)

sparse.baps <- sparse.base
sparse.baps$snp.matrix <- sparse.baps$snp.matrix[, colnames(sparse.baps$snp.matrix) %in% raxml_tree$tip.label]
sparse.baps <- optimise_prior(sparse.baps, type = "baps", hc.method = "ward", grid.interval = c(5e-9, 10), n.cores = 1)
d.baps <- snp_dist(sparse.baps)
d.baps <- as.dist(d.baps / max(d.baps))
h.baps <- hclust(d.baps, method = "ward.D2")
partition.baps <- best_baps_partition(sparse.baps, h.baps)
partition.baps <- as.data.frame(partition.baps) %>% rownames_to_column()
colnames(partition.baps) <- c('tip', 'group.baps')
partition.baps$group.baps <- as.factor(partition.baps$group.baps)

ggtree(raxml_tree) %<+% partition.baps + geom_tippoint(aes(color = group.baps), size = 2) + 
  paletteer::scale_color_paletteer_d("ggsci::category20_d3")
  # paletteer::scale_color_paletteer_d("pals::polychrome")
  # paletteer::scale_color_paletteer_d("pals::alphabet")

# 查看配色数量
# paletteer::paletteer_d("pals::alphabet") %>% length()

boot.baps   <- boot_fast_baps(sparse.baps, n.cores = 8, quiet = FALSE)
dendro.baps <- stats::as.dendrogram(h.baps)
gplots::heatmap.2(boot.baps, dendro.baps, dendro.baps, tracecol = NA)


sparse.opt <- sparse.base
sparse.opt$snp.matrix <- sparse.opt$snp.matrix[, colnames(sparse.opt$snp.matrix) %in% raxml_tree$tip.label]
sparse.opt <- optimise_prior(sparse.opt, type = "optimise.baps", hc.method = "ward", grid.interval = c(5e-9, 10), n.cores = 1)
d.opt <- snp_dist(sparse.opt)
d.opt <- as.dist(d.opt / max(d.opt))
h.opt <- hclust(d.opt, method = "ward.D2")
partition.opt <- best_baps_partition(sparse.opt, h.opt)
partition.opt <- as.data.frame(partition.opt) %>% rownames_to_column()
colnames(partition.opt) <- c('tip', 'group.opt')
partition.opt$group.opt <- as.factor(partition.opt$group.opt)
rm(d, h)

ggtree(raxml_tree) %<+% partition.opt + geom_tippoint(aes(color = group.opt), size = 2) + 
  paletteer::scale_color_paletteer_d("ggsci::category20_d3")

boot.opt   <- boot_fast_baps(sparse.opt, n.cores = 8, quiet = FALSE)
dendro.opt <- stats::as.dendrogram(h.opt)
gplots::heatmap.2(boot.opt, dendro.opt, dendro.opt, tracecol = NA)


sparse.tree <- sparse.base
sparse.tree$snp.matrix <- sparse.tree$snp.matrix[, colnames(sparse.tree$snp.matrix) %in% raxml_tree$tip.label]
sparse.tree <- optimise_prior(sparse.tree, type = "optimise.baps", hc.method = "ward", grid.interval = c(5e-9, 10), n.cores = 1)
partition.tree <- best_baps_partition(sparse.tree, raxml_tree)
partition.tree <- as.data.frame(partition.tree) %>% rownames_to_column()
colnames(partition.tree) <- c('tip', 'group.tree')
partition.tree$group.tree <- as.factor(partition.tree$group.tree)

ggtree(raxml_tree) %<+% partition.tree + geom_tippoint(aes(color = group.tree), size = 2) + 
  paletteer::scale_color_paletteer_d("ggsci::category20_d3")

boot.tree   <- boot_fast_baps(sparse.tree, n.cores = 8, quiet = FALSE)
dendro.tree <- phylogram::as.dendrogram(ape::compute.brlen(raxml_tree, method = "Grafen"))
gplots::heatmap.2(boot.tree, dendro.tree, dendro.tree, tracecol = NA)



p <- ggtree(raxml_tree)
p <- facet_plot(p, panel = "fastBAPS unoptimized\nInitial clusters based on SNPs", data = partition.baps, geom = geom_tile, aes(x = as.integer(group.baps)), color = '#7f9faf')
p <- facet_plot(p, panel = "fastBAPS optimized\nInitial clusters based on SNPs", data = partition.opt, geom = geom_tile, aes(x = as.integer(group.opt)), color = '#cf929e')
p <- facet_plot(p, panel = "fastBAPS optimized\nInitial clusters based on RAxML tree", data = partition.tree, geom = geom_tile, aes(x = as.integer(group.tree)), color = '#06436f')
p

save.image('./r_image/clustering_analysis_fastBAPS.RData')
