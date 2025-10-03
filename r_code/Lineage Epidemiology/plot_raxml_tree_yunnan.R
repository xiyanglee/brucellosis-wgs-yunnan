library(Biostrings) # 读取 fasta 文件
library(tidyverse)
library(ape)
library(phangorn)
library(ggtree)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



# 读取数据
fasta_file <- readDNAStringSet("./data/snippy_contig_gubbins/clean.core.aln")

# 设置外类群
# "GCA_000369945.1"
outgroup_tip <- "GCA_000022625.1" 
all_tips <- union(metadata.sample$assembly_accession, outgroup_tip)
matched_tips <- intersect(all_tips, names(fasta_file))

# 提取序列
aligned_seqs <- fasta_file[matched_tips]

# 转换为 phangorn 可用格式：phyDat
aligned_seqs_phydat <- phyDat(as.matrix(aligned_seqs), type = "DNA")

# 构建初始 NJ 树作为 ML 起点
dm <- dist.ml(aligned_seqs_phydat)
tree_start <- NJ(dm)

# 建立 ML 模型
fit <- pml(tree_start, data = aligned_seqs_phydat)

# 优化树（可加 bootstrap、模型选择）
# “JC”, “F81”, “K80”, “HKY”, “TrNe”, “TrN”, “TPM1”, “K81”, “TPM1u”, “TPM2”, “TPM2u”, “TPM3”, “TPM3u”, 
# “TIM1e”, “TIM1”, “TIM2e”, “TIM2”, “TIM3e”, “TIM3”, “TVMe”, “TVM”, “SYM”, “GTR”
fit_opt <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic")

# 提取优化后的树
tree_ml <- fit_opt$tree

# 用外类群定根
tree_ml_rooted <- root(tree_ml, outgroup = outgroup_tip, resolve.root = TRUE)
# 裁剪外类群后绘图
tree_pruned <- drop.tip(tree_ml_rooted, outgroup_tip)

# 先画一次树并展示所有节点编号（含内部节点）
ggtree(tree_pruned) +
  geom_text(aes(label = node), size = 3, color = "blue") +
  labs(title = "Tree with Node IDs (for rooting reference)") +
  theme_tree2()

tree_pruned <- root(tree_pruned, node = 176, resolve.root = TRUE)

# 关联分组信息
tree_data <- data.frame(tip.label = tree_pruned$tip.label) %>%
  left_join(., metadata.sample, by = c("tip.label" = "assembly_accession"))

# ggtree(tree_pruned) %<+% tree_data + 
ggtree(ape::compute.brlen(tree_pruned, method = "Grafen"), size = 0.3, layout = "fan", open.angle = 180) %<+% tree_data + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = year(collection_date)),
    width = 0.05,
    offset = 0.1
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = age),
    width = 0.05,
    offset = 0.1
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = gender),
    width = 0.05,
    offset = 0.1
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = patient_type),
    width = 0.05,
    offset = 0.1
  ) + 
  geom_hilight(node = 175, alpha = 0.3, fill = '#7EBB8F') +
  geom_hilight(node = 176, alpha = 0.3, fill = '#9DD0CF') +
  geom_hilight(node = 193, alpha = 0.3, fill = '#886A9E') +
  geom_cladelabel(node = 175, label = "1.1", align = T, color = '#5e8c6b', linewidth = 2) +
  geom_cladelabel(node = 176, label = "1.2", align = T, color = '#7da6a5', linewidth = 2) + 
  geom_cladelabel(node = 193, label = "1.3", align = T, color = '#6e5680', linewidth = 2)

ggtree(
  ape::compute.brlen(tree_pruned, method = "Grafen"),
  size = 0.75, layout = "fan", open.angle = 180
) %<+% 
  (
    tree_data %>% 
      mutate(
        patient_type = dplyr::recode(
          patient_type,
          '其他省' = 'out-of-city',
          '本省其它地市' = 'out-of-city',
          '本县区' = 'in-city',
          '本市其它县区' = 'in-city'
        ),
        patient_type = factor(patient_type, levels = c('in-city', 'out-of-city')) 
      )
  ) +
  # 1) collection year (continuous) 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = year(collection_date)),
    width = 0.075, offset = 0.09
  ) +
  scale_fill_gradient(
    low  = "#aed0ee", 
    high = "#003d74", 
    name = "Year"
  ) +
  ggnewscale::new_scale_fill() +
  
  # 2) age (continuous) 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = age),
    width = 0.075, offset = 0.09
  ) +
  scale_fill_gradient(
    low  = "#ffeaa7", 
    high = "#d23918", 
    name = "Age"
  ) +
  ggnewscale::new_scale_fill() +
  
  # 3) gender (discrete) 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = gender),
    width = 0.075, offset = 0.09
  ) +
  scale_fill_manual(
    values = c(
      "男" = "#dcc7e1", 
      "女" = "#7d5284" 
    ),
    name = "Gender"
  ) +
  ggnewscale::new_scale_fill() +
  
  # 4) patient type (discrete)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip.label, fill = patient_type),
    width = 0.075, offset = 0.09
  ) +
  scale_fill_manual(
    values = c(
      "out-of-city" = "#778ca3",
      "in-city"     = "#e0e0d0"
    ),
    name = "Patient Type"
  ) +
  
  theme(legend.position = "bottom") +
  
  # highlights and clade labels
  geom_hilight(node = 175, alpha = 0.4, fill = "#66aa88") +
  geom_hilight(node = 193, alpha = 0.4, fill = "#88cccc") +
  geom_hilight(node = 176, alpha = 0.4, fill = "#775599") +
  geom_cladelabel(node = 175, label = "1.1", align = TRUE,  color = "#3e6b47", linewidth = 2) +
  geom_cladelabel(node = 193, label = "1.2", align = TRUE,  color = "#4d8b8b", linewidth = 2) +
  geom_cladelabel(node = 176, label = "1.3", align = TRUE,  color = "#5e3b6a", linewidth = 2)

# ggsave('./plot/tree/raxml_tree_snippy_yunnan_v2_type2.pdf', width = 8, height = 8)



