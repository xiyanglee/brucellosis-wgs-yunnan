# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("gtonkinhill/fastbaps", build_vignettes = TRUE))
# withr::with_libpaths("/public/r_share_library/4.1", remotes::install_github('YuLab-SMU/ggtree'))
# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("xiangpin/ggtreeExtra"))
library(fastbaps)
library(ape)
library(Matrix)
library(tidyverse)
library(ggtree)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

raxml_tree <- ape::read.tree('./data/core_genes_removed_gap/RAxML_bestTree.raxml_out')


# baps.hc <- fast_baps(fasta_sparse_data)
# best.partition <- best_baps_partition(fasta_sparse_data, as.phylo(baps.hc))

fasta2tree_name_comparison_table <- read.table('./data/core_genes_removed_gap/phylip_sample_name_genebank.tsv', sep = '\t', header = FALSE, 
                                               col.names = c('fasta', 'tree'), colClasses = 'character')
raxml_tree$tip.label <- plyr::mapvalues(raxml_tree$tip.label, from = fasta2tree_name_comparison_table$tree, to = fasta2tree_name_comparison_table$fasta)


plot(raxml_tree, show.node.label = T)

# "phylo" 对象中，叶节点编号为 1 到 Ntip(tree)，内部节点编号为 Ntip(tree) + 1 到 Nnode(tree) + Ntip(tree)
tip_depths <- node.depth.edgelength(raxml_tree)[1:Ntip(raxml_tree)]
tip_depths_df <- data.frame(tip = raxml_tree$tip.label, depth = tip_depths)
tip_depths_df$z_score <- (tip_depths_df$depth - mean(tip_depths_df$depth)) / sd(tip_depths_df$depth)

drop_tip <- tip_depths_df$tip[tip_depths_df$z_score > 1.7]
raxml_tree <- drop.tip(raxml_tree, drop_tip)



fasta_sparse_data <- import_fasta_sparse_nt('./data/core_genes_removed_gap/core_gene_snp_filtered.aln')

keep_tip <- !(colnames(fasta_sparse_data$snp.matrix) %in% drop_tip) & (colnames(fasta_sparse_data$snp.matrix) %in% raxml_tree$tip.label)
fasta_sparse_data$snp.matrix <- fasta_sparse_data$snp.matrix[, keep_tip]

fasta_sparse_data <- optimise_prior(fasta_sparse_data, type = "optimise.symmetric", hc.method = "ward", grid.interval = c(5e-9, 1))
# [1] "Optimised hyperparameter: 0.007"


raxml_hclust <- phytools::midpoint.root(raxml_tree)
best.partition <- best_baps_partition(fasta_sparse_data, raxml_hclust)


group_data <- as.data.frame(best.partition) %>% 
  rownames_to_column()
colnames(group_data) <- c('tip', 'group')
group_data$group <- as.factor(group_data$group)

# ggtree(raxml_tree, layout = 'equal_angle') %<+% metadata + geom_tippoint(aes(color = group), size = 2)
ggtree(raxml_tree) %<+% group_data + geom_tippoint(aes(color = group), size = 2) + 
  ggsci::scale_color_bmj()
# ggtree(as.phylo(baps.hc)) %<+% metadata + geom_tippoint(aes(color = group), size = 2)



load('./r_image/metadata.RData')

tip_depths_df <- dplyr::filter(tip_depths_df, !tip %in% drop_tip)

repeat_data <- read.csv('./data/core_genes_removed_gap/core_gene_alignment_repeat.csv')
repeat_data <- left_join(repeat_data, tip_depths_df, by = c('target_seqname' = 'tip'))
repeat_data$target_seqname <- NULL
colnames(repeat_data)[1] <- 'tip'
tip_depths_all <- rbind(repeat_data, tip_depths_df)

repeat_data <- read.csv('./data/core_genes_removed_gap/core_gene_alignment_repeat.csv')
repeat_data <- left_join(repeat_data, group_data, by = c('target_seqname' = 'tip'))
repeat_data$target_seqname <- NULL
colnames(repeat_data)[1] <- 'tip'
group_data_all <- rbind(repeat_data, group_data)

metadata <- left_join(left_join(tip_depths_all, group_data_all, by = 'tip'), metadata, by = c('tip' = 'assembly_accession'))

ggplot(filter(metadata, group == 2), aes(x = collection_year, y = depth, color = area)) + 
  geom_point()


metadata <- metadata %>% 
  mutate(
    label_yunnan = case_when(
      continent == 'Yunnan' ~ 'Yunnan', 
      .default = NA), 
    continent = case_when(
      continent == 'Yunnan' ~ 'Asia', 
      .default = continent
    )
  )

library(ggtreeExtra)
ggtree(ape::compute.brlen(raxml_tree, method = "Grafen"), size = 0.3, layout = "fan", open.angle = 180) %<+% 
  metadata + 
  # geom_tippoint(aes(color = group), size = 2) +
  # geom_tiplab(size = 3, hjust = 0.2) +
  # theme_minimal() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = label_yunnan),
    width = 0.1,
    offset = 0.1
  ) + 
  scale_fill_manual(values = c("Yunnan" = "gold"), na.value = "transparent") + # 将 NA 值设置为透明 
  ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = group),
    width = 0.1,
    offset = 0.1
  ) + 
  scale_fill_manual(values = c("1" = "#BEB8DC", "2" = "#F7E1ED", "3" = "#CFEAF1"), na.value = "transparent") + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = continent),
    width = 0.1,
    offset = 0.1
  ) + 
  scale_fill_manual(values = c("Asia" = "#8ECFC9", "Europe" = "#FFBE7A", "Africa" = "#82B0D2", "America" = "#FA7F6F"), 
                    na.value = "grey95") + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip, fill = collection_year),
    width = 0.1,
    offset = 0.1
  ) + 
  # scale_fill_gradientn(colors = c("#2878B5", "#E7EFFA", "#F3D266", "#C82423")) 
  scale_fill_gradientn(colors = c("#496C88", "#E7EFFA", "#F8F3F9", "#FF8884")) 
  # ggsci::scale_color_bmj()

ggsave('./plot/tree/raxml_tree.pdf', width = 8, height = 8)

