library(tidyverse)
library(Biostrings)
library(ape)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata.RData')

metadata_filter <- metadata[!is.na(metadata$collection_date) & !is.na(metadata$area), ]

# 表格 metadata_filter，其中有三列 "collection_year"  "collection_month" "collection_day"，
# 其中 "collection_month" "collection_day" 可能为 NA，
# 
# 编写函数生成一个新列 "collection_date"，使用浮点数表示，
# 如 2020..002739726 表示 2020 年 1 月 1 日；
# 计算方法为，若"collection_year"  "collection_month" "collection_day"都不为 NA，则 collection_date 为对应日期为当年第多少天再除以 365，
# 若 "collection_day" 为 NA，则 collection_date =  collection_year + collection_month /12，
# 若 "collection_month" "collection_day" 都为空，则 collection_date = collection_year；
# 
# 同时生成一列 “uncertainty”，若 "collection_month" "collection_day" 均为 NA，则 uncertainty = 1，
# 若 "collection_day" 为空，则 uncertainty = 1/12，其他情况 uncertainty = 0；

metadata_filter <- metadata_filter %>%
  mutate(
    # 计算 collection_date
    collection_date = case_when(
      is.na(collection_month) & is.na(collection_day) ~ as.numeric(collection_year),
      is.na(collection_day) ~ as.numeric(collection_year) + collection_month / 12,
      !is.na(collection_month) & !is.na(collection_day) ~ {
        day_of_year <- as.numeric(format(as.Date(
          paste(collection_year, collection_month, collection_day, sep = "-"),
          format = "%Y-%m-%d"
        ), "%j"))
        collection_year + (day_of_year - 1) / 365
      }
    ),
    # 计算 uncertainty
    uncertainty = case_when(
      is.na(collection_month) & is.na(collection_day) ~ 1,
      is.na(collection_day) ~ 1 / 12,
      .default = 0
    )
  )

metadata_filter <- metadata_filter %>% dplyr::select(assembly_accession, collection_date, geographic_country, continent, area, uncertainty)
write.csv(metadata_filter, './data/beast/metadata.csv', row.names = FALSE)

metadata_filter <- filter(metadata_filter, ! assembly_accession %in% c('GCA_000250835.2', 'GCA_039649315.1', 'GCA_009823695.1', 'GCA_000250775.2'))
write.csv(metadata_filter, './data/beast_gubbins_contig/metadata.csv', row.names = FALSE)



# alignment <- read.dna('./data/core_genes_removed_gap/core_gene_snp_filtered.aln', format = 'fasta')
alignment <- read.dna("./data/snippy_contig_gubbins/clean.core.aln", format = "fasta")
filtered_alignment <- alignment[metadata_filter$assembly_accession, ]

# 保存为 NEXUS 文件
# write.nexus.data(filtered_alignment, file = "./data/beast/core_gene_snp_filtered.nex", format = "DNA")
write.nexus.data(filtered_alignment, file = "./data/beast_gubbins_contig/clean_core_filtered.nex", format = "DNA")
