library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')




normalize_datamat <- function(datmat) {
  colnames(datmat) <- datmat[1, ]
  datmat <- datmat[-1, ]
  
  rownames_datamat <- datmat$Gene
  datmat$Gene <- NULL
  
  datmat <- datmat %>% reframe(across(everything(), as.integer))
  rownames(datmat) <- rownames_datamat
  
  return(datmat)
}

pangenes.panaroo  <- read.csv('./data/panaroo_result/gene_presence_absence_roary.csv', header = TRUE, stringsAsFactors = FALSE)
gene_data.panaroo <- read.csv('./data/panaroo_result/gene_data.csv', header = TRUE, stringsAsFactors = FALSE)
gene_data.panaroo <- dplyr::select(gene_data.panaroo, -c(prot_sequence, dna_sequence))

datmat.panaroo <- read.table('./data/panaroo_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat.panaroo <- normalize_datamat(datmat.panaroo)
# datmat.roary <- read.table('./data/roary_result/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
# datmat.roary <- normalize_datamat(datmat.roary)

head(pangenes.panaroo[, c("Gene", "Non.unique.Gene.name")])
#            Gene Non.unique.Gene.name
# 1        flgG_1               flgG_1
# 2 nanR_1~~~gntR         ;nanR_1;gntR
# 3    group_1428                     
# 4    group_1427                     
# 5 xerC_1~~~xerC          xerC_1;xerC
# 6    group_1425                     

gene_family_map <- pangenes.panaroo %>%
  filter(Non.unique.Gene.name != "") %>%   # 去除无 gene name 的行
  select(family_name = Gene, gene_name = Non.unique.Gene.name) %>%  # 重命名列
  separate_rows(gene_name, sep = ";") %>%  # 按 ; 拆分为多行
  filter(gene_name != "")                  # 去除空 gene_name

gene_family_map <- gene_family_map %>%
  group_by(family_name) %>%
  filter(n_distinct(gene_name) > 1) %>%
  ungroup()

family_gene_names <- unique(gene_family_map$gene_name)

sum(!family_gene_names %in% gene_data.panaroo$gene_name)
sum(gene_data.panaroo$gene_name == "")

sum(str_starts(gene_data.panaroo$gene_name, "group_"))
sum(str_starts(family_gene_names, "group_"))


gene_data.panaroo <- gene_data.panaroo[gene_data.panaroo$gene_name %in% family_gene_names, ]
gene_data.panaroo <- left_join(gene_data.panaroo, gene_family_map, by = "gene_name")
gene_data.panaroo <- gene_data.panaroo[!is.na(gene_data.panaroo$family_name), ]

datmat.gene_family <- as.data.frame(table(gene_data.panaroo$gff_file, gene_data.panaroo$family_name))
datmat.gene_family <- datmat.gene_family %>%
  pivot_wider(
    names_from = Var2,
    values_from = Freq,
    values_fill = 0
  ) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "Var1")

family_genes_count <- as.data.frame(rowMeans(datmat.gene_family)) %>% rownames_to_column("assembly_accession")
colnames(family_genes_count)[2] <- "gene_counts"






beast_tree <- treeio::read.beast('./data/beast_gubbins_contig/run_1/beast_GTRIG_clean_core_filtered.tree')

MAX_DATE <- 2023.58333333333
mrca.height <- beast_tree[, c("node", "height")]
mrca.height <- mrca.height %>% mutate(height = as.numeric(height))

# ---------------------------
# 生成两两样本对，并获取样本对的 MRCA 编号
# 
sample_ids <- intersect(metadata.phylo$assembly_accession, beast_tree@phylo$tip.label)
sample_pairs <- combn(sample_ids, 2)

# 获取 MRCA 编号
mrcas <- vapply(
  seq_len(ncol(sample_pairs)),
  function(i) getMRCA(beast_tree@phylo, sample_pairs[, i]),
  FUN.VALUE = integer(1)
)

mrca.nodes <- data.frame(
  sample_1 = sample_pairs[1, ],
  sample_2 = sample_pairs[2, ],
  mrca_node = mrcas
)

mrca.nodes <- left_join(mrca.nodes, mrca.height, by = c('mrca_node' = 'node'))

rm(sample_pairs, mrcas)

mrca.nodes <- mrca.nodes %>%
  left_join(., family_genes_count, by = c("sample_1" = "assembly_accession")) %>%
  dplyr::rename(gene_count_1 = gene_counts) %>%
  left_join(., family_genes_count, by = c("sample_2" = "assembly_accession")) %>%
  dplyr::rename(gene_count_2 = gene_counts) %>%
  mutate(diff_gene_counts = abs(gene_count_1 - gene_count_2))

plot_df <- mrca.nodes[, c("height", "diff_gene_counts")]
plot_df$diff_gene_counts <- plot_df$diff_gene_counts * 1161
# plot_df$log10_height <- log10(plot_df$height)
# plot_df$sqrt_diff_gene_counts <- sqrt(plot_df$diff_gene_counts)

# 拟合百分位回归模型
fit_q025 <- quantreg::rq(diff_gene_counts ~ height - 1, tau = 0.025, data = plot_df)
fit_q50  <- quantreg::rq(diff_gene_counts ~ height - 1, tau = 0.5,   data = plot_df)
fit_q975 <- quantreg::rq(diff_gene_counts ~ height - 1, tau = 0.975, data = plot_df)
# fit_q025 <- quantreg::rq(sqrt_diff_gene_counts ~ log10_height - 1, tau = 0.025, data = plot_df)
# fit_q50  <- quantreg::rq(sqrt_diff_gene_counts ~ log10_height - 1, tau = 0.5,   data = plot_df)
# fit_q975 <- quantreg::rq(sqrt_diff_gene_counts ~ log10_height - 1, tau = 0.975, data = plot_df)

# fit_q025$coefficients
# height 
# 0.0008049982 
# fit_q50$coefficients
# height 
# 0.006171653 
# fit_q975$coefficients
# height 
# 0.08667203

# 为每个 height 值预测相应的 diff_gene_counts 百分位回归值
plot_df <- plot_df %>%
  mutate(
    Q2.5  = predict(fit_q025, newdata = .),
    Q50   = predict(fit_q50,  newdata = .),
    Q97.5 = predict(fit_q975, newdata = .)
    # Q2.5  = Q2.5 ** 2, 
    # Q50   = Q50 ** 2, 
    # Q97.5 = Q97.5 ** 2
  )

library(ggrastr)

ggplot(plot_df, aes(x = height, y = diff_gene_counts)) + 
  ## -> 1. 正常图像
  # geom_point(color = "#b2bfc3", alpha = 0.05, size = 3) + 
  ## -> 2. rasterise, 减小文件体积
  ## -> 3. rasterise 方法在浅色时会偏色，可以使用正常导出之后在 AI 中 "拼合透明度" 减小文件体积
  rasterise(geom_point(color = "#a8afc3", alpha = 0.05, size = 3), dpi = 600) + 
  geom_line(aes(y = Q2.5),  color = "black", linetype = "dashed", linewidth = 0.75) +
  geom_line(aes(y = Q50),   color = "black", linewidth = 1) +
  geom_line(aes(y = Q97.5), color = "black", linetype = "dashed", linewidth = 0.75) +
  scale_x_log10(name = "Evolutionary time (years)") + 
  # scale_y_log10() + 
  scale_y_sqrt(breaks = c(0, 10, 100, 200, 300, 400, 500), 
               name = "Paralogs copy number variation",
               sec.axis = sec_axis(
                 transform = ~ . / 1161,
                 breaks = c(0, 10, 100, 200, 300, 400, 500) / 1161,
                 labels = scales::number_format(accuracy = 0.01),
                 name = "Mean per-gene copy number variation"
               )) +
  annotation_logticks(sides = "b") +
  geom_segment(data = data.frame(y = c(seq(0, 10, 1), seq(20, 100, 10), seq(200, 1000, 100))),
               aes(y = y, yend = y, x = 1, xend = 1.1),
               inherit.aes = FALSE) +
  geom_segment(data = data.frame(y = c(5, 50, 500)),
               aes(y = y, yend = y, x = 1, xend = 1.2),
               inherit.aes = FALSE) +
  geom_segment(data = data.frame(y = c(10, 100, 1000)),
               aes(y = y, yend = y, x = 1, xend = 1.3),
               inherit.aes = FALSE) +
  coord_cartesian(xlim = c(1, NA), ylim = c(0, 550), expand = 0) +
  theme_bw(base_size = 14)

# ggsave('./plot/pan_genome_analysis/homologous_copy_number_variation.pdf', width = 4.3, height = 3.5)
# ggsave('./plot/pan_genome_analysis/homologous_copy_number_variation_rasterise.pdf', width = 4.3, height = 3.5)

write.csv(plot_df, "./r_plot_data/figure_1e~plot_df.csv", row.names = FALSE)
