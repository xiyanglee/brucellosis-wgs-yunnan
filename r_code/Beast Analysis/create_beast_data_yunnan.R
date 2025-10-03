library(ape)

load('./r_image/metadata.RData')
load('./r_image/cluster_tree.RData')

alignment <- read.dna("./data/snippy_contig_gubbins/clean.core.aln", format = "fasta")
filtered_alignment <- alignment[metadata.sample$assembly_accession, ]
rownames(filtered_alignment) <- paste(metadata.sample$assembly_accession, metadata.sample$collection_date, sep = '|')

# 保存为 NEXUS 文件
# write.nexus.data(filtered_alignment, file = "./data/beast/core_gene_snp_filtered.nex", format = "DNA")
write.nexus.data(filtered_alignment, file = "./data/beast_gubbins_contig_yunnan/clean_core_filtered.nex", format = "DNA")

metadata.sample$traits <- paste(metadata.sample$assembly_accession, metadata.sample$collection_date, sep = '|')

seq_loc <- read.csv('../../metadata/Yunnan_sequenced_patient_location.csv')
seq_loc$input_address <- NULL
colnames(seq_loc) <- c('assembly_accession', 'long', 'lat')
seq_loc$assembly_accession <- as.character(seq_loc$assembly_accession)

metadata.sample <- left_join(seq_loc, metadata.sample, by = 'assembly_accession')
metadata.sample <- metadata.sample %>% dplyr::select(traits, lat, long)

write.table(metadata.sample, './data/beast_gubbins_contig_yunnan/coordinates.txt', sep = '\t', row.names = FALSE, quote = FALSE)



load('./r_image/metadata.RData')
load('./r_image/cluster_tree.RData')
metadata.sample <- as.data.frame(metadata.sample)
rownames(metadata.sample) <- metadata.sample$assembly_accession

alignment <- read.dna("./data/snippy_contig_gubbins/clean.core.aln", format = "fasta")

sample_name_1.1 <- group_data_2$tip[group_data_2$group_2 == '1.1']
sample_name_1.1 <- sample_name_1.1[!grepl("^GCA_", sample_name_1.1)]
filtered_alignment_1.1 <- alignment[sample_name_1.1, ]
rownames(filtered_alignment_1.1) <- paste(rownames(filtered_alignment_1.1), metadata.sample[rownames(filtered_alignment_1.1), 'collection_date'], sep = '|')
write.nexus.data(filtered_alignment_1.1, file = "./data/beast_gubbins_contig_yunnan_group/1.1/clean_core_filtered.nex", format = "DNA")

sample_name_1.2 <- group_data_2$tip[group_data_2$group_2 == '1.2']
sample_name_1.2 <- sample_name_1.2[!grepl("^GCA_", sample_name_1.2)]
filtered_alignment_1.2 <- alignment[sample_name_1.2, ]
rownames(filtered_alignment_1.2) <- paste(rownames(filtered_alignment_1.2), metadata.sample[rownames(filtered_alignment_1.2), 'collection_date'], sep = '|')
write.nexus.data(filtered_alignment_1.2, file = "./data/beast_gubbins_contig_yunnan_group/1.2/clean_core_filtered.nex", format = "DNA")

sample_name_1.3 <- group_data_2$tip[group_data_2$group_2 == '1.3']
sample_name_1.3 <- sample_name_1.3[!grepl("^GCA_", sample_name_1.3)]
filtered_alignment_1.3 <- alignment[sample_name_1.3, ]
rownames(filtered_alignment_1.3) <- paste(rownames(filtered_alignment_1.3), metadata.sample[rownames(filtered_alignment_1.3), 'collection_date'], sep = '|')
write.nexus.data(filtered_alignment_1.3, file = "./data/beast_gubbins_contig_yunnan_group/1.3/clean_core_filtered.nex", format = "DNA")

metadata.sample$traits <- paste(metadata.sample$assembly_accession, metadata.sample$collection_date, sep = '|')
seq_loc <- read.csv('../../metadata/Yunnan_sequenced_patient_location.csv')
seq_loc$input_address <- NULL
colnames(seq_loc) <- c('assembly_accession', 'long', 'lat')
seq_loc$assembly_accession <- as.character(seq_loc$assembly_accession)
metadata.sample <- left_join(seq_loc, metadata.sample, by = 'assembly_accession')
metadata.sample <- left_join(metadata.sample, group_data_2, by = c('assembly_accession' = 'tip'))

metadata.sample.1.1 <- dplyr::filter(metadata.sample, group_2 == '1.1')
metadata.sample.1.1 <- metadata.sample.1.1 %>% dplyr::select(traits, lat, long)
write.table(metadata.sample.1.1, './data/beast_gubbins_contig_yunnan_group/1.1/coordinates.txt', sep = '\t', row.names = FALSE, quote = FALSE)

metadata.sample.1.2 <- dplyr::filter(metadata.sample, group_2 == '1.2')
metadata.sample.1.2 <- metadata.sample.1.2 %>% dplyr::select(traits, lat, long)
write.table(metadata.sample.1.2, './data/beast_gubbins_contig_yunnan_group/1.2/coordinates.txt', sep = '\t', row.names = FALSE, quote = FALSE)

metadata.sample.1.3 <- dplyr::filter(metadata.sample, group_2 == '1.3')
metadata.sample.1.3 <- metadata.sample.1.3 %>% dplyr::select(traits, lat, long)
write.table(metadata.sample.1.3, './data/beast_gubbins_contig_yunnan_group/1.3/coordinates.txt', sep = '\t', row.names = FALSE, quote = FALSE)
