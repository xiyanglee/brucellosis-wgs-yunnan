library(seqinr)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')



raxml_tree <- ape::read.tree('./data/snippy_contig_gubbins/gubbins.final_tree.tre')
plot(raxml_tree, show.node.label = T)

drop_tip <- c('GCA_000369945.1', 'GCA_000250835.2', 'GCA_039649315.1', 'GCA_009823695.1', 'GCA_000250775.2')

raxml_tree <- drop.tip(raxml_tree, drop_tip)
plot(raxml_tree, show.node.label = T)



alignment <- read.alignment(file = "./data/snippy_contig_gubbins/clean.core.aln", format = "fasta")
alignment$seq <- lapply(alignment$seq, toupper)

# 去除指定名称的序列
remove_names <- c('GCA_000250835.2', 'GCA_009823695.1', 'GCA_000250775.2')
keep_indices <- !(alignment$nam %in% remove_names)
filtered_alignment <- list(
  nb = sum(keep_indices),
  nam = alignment$nam[keep_indices],
  seq = alignment$seq[keep_indices]
)

# 写入过滤后的文件
write.fasta(sequences = filtered_alignment$seq, 
            names = filtered_alignment$nam, 
            file.out = "./data/raxml_gubbins_contig/filtered_branch.core.aln")



# alignment <- read.alignment(file = "./data/snippy_contig_gubbins/clean.core.aln", format = "fasta")
# alignment$seq <- lapply(alignment$seq, toupper)
# 
# # 去除指定名称的序列
# remove_names <- c('GCA_000250835.2', 'GCA_039649315.1', 'GCA_009823695.1', 'GCA_000250775.2')
# keep_indices <- !(alignment$nam %in% remove_names)
# filtered_alignment <- list(
#   nb = sum(keep_indices),
#   nam = alignment$nam[keep_indices],
#   seq = alignment$seq[keep_indices]
# )
# 
# # 写入过滤后的文件
# write.fasta(sequences = filtered_alignment$seq, 
#             names = filtered_alignment$nam, 
#             file.out = "./data/beast_gubbins_contig/filtered_outgroup_branch.core.aln")
