# withr::with_libpaths("/public/r_share_library/4.1", devtools::install_github("xavierdidelot/BactDating"))
library(tidyverse)
library(BactDating)
library(ape)
library(coda)

setwd('/Users/xiyangli/Lab/Project/Brucella_WGS_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



metadata.phylo <- metadata.phylo %>%
  mutate(
    # 计算 collection_date
    collection_date_float = case_when(
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
      is.na(collection_date_float) ~ NA, 
      is.na(collection_month) & is.na(collection_day) ~ 1,
      is.na(collection_day) ~ 1 / 12,
      .default = 0
    )
  )

aln <- read.dna("./data/raxml_gubbins_contig/filtered_branch.core.aln.reduced", format = "sequential", as.character = TRUE)
# aln <- read.dna("./data/raxml_gubbins_contig/filtered_branch.core.aln", format = "fasta", as.character = TRUE)
aln_length <- ncol(as.matrix(aln))
rm(aln)

raxml_tree$edge.length <- raxml_tree$edge.length * aln_length

# 创建时间上下限矩阵（可能带 NA）
d_matrix <- with(metadata.phylo, 
                 cbind(collection_date_float, 
                       collection_date_float + uncertainty))
# 设置行名为样本 ID（必须与树 tip.label 匹配）
rownames(d_matrix) <- metadata.phylo$assembly_accession

# 4. 使用 BactDating（指定 relaxed clock 模型、gamma prior）
dated_tree <- bactdate(tree = raxml_tree,
                       date = d_matrix,
                       model = "mixedgamma",    # 混合 Gamma 模型
                       nbIts = 1000000,         # MCMC 迭代数
                       updateRoot = FALSE,      # 是否估计根节点时间
                       showProgress = TRUE)

plot(dated_tree, 'treeCI', show.tip.label = F)
plot(dated_tree, 'treeRoot', show.tip.label = F)
plot(dated_tree, 'trace')

plot(extractSample(dated_tree, 100)[[20]])
axisPhylo(1,backward = F)




effectiveSize(mcmc)
# save.image('./r_data/BactDating.RData')

