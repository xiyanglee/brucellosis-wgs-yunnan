# BiocManager::install("clusterProfiler", lib = "/public/r_share_library/4.1")
# BiocManager::install("KEGGREST", lib = "/public/r_share_library/4.1")
library(clusterProfiler)
library(KEGGREST)



# QUERY2GENE <- readLines('./data/roary_result/pan_genome_reference.fa') %>% 
#   grep("^>", ., value = TRUE) %>% 
#   sub("^>", "", .) %>% 
#   strsplit(" ")
# 
# QUERY2GENE <- data.frame(
#   QUERY = sapply(QUERY2GENE, `[`, 1), 
#   GENE  = sapply(QUERY2GENE, `[`, 2)
# )



readDelimWithDoubleChar <- function(x, ...) {
  txt <- readLines(x)
  txt <- sub('^#([^#]*)', '\\1', txt)
  read.delim(text = txt, comment.char = "#", ...)
}

# emapper_annot <- readDelimWithDoubleChar('./data/panaroo_emapper_diamond/pan_genome_annotation.emapper.annotations', sep = '\t')
emapper_annot <- readDelimWithDoubleChar('./data/panaroo_result_strict_pan/pan_genome_annotation.emapper.annotations', sep = '\t')
emapper_annot[emapper_annot == "-"] <- NA 
emapper_annot <- as.data.frame(emapper_annot)
# emapper_annot <- left_join(emapper_annot, QUERY2GENE, by = c('query' = 'QUERY'))



## build KO orthology BY emapper predicted
## --> db.KO.emapper

emapper.Protein2KO <- emapper_annot %>% 
  dplyr::select(query, KEGG_ko) %>% 
  tidyr::drop_na(KEGG_ko)

emapper.TERM2GENE <- NULL
for (row in 1:nrow(emapper.Protein2KO)) {
  TERM <- str_split(emapper.Protein2KO[row, "KEGG_ko"], ",", simplify = FALSE)[[1]]  
  GENE <- emapper.Protein2KO[row, "query"][[1]]
  tmp <- data.frame(TERM = TERM, GENE = GENE)
  emapper.TERM2GENE <- rbind(emapper.TERM2GENE, tmp)
}; rm(row, GENE, tmp, TERM)
emapper.TERM2GENE <- unique(emapper.TERM2GENE)
# emapper.TERM2NAME <- ko2name(unique(emapper.TERM2GENE$TERM))
# ko2name 函数结果不正确，根据 API 重新获取 KO 和 name 的对应关系
KEGG.ko_list <- KEGGREST:::.getUrl('https://rest.kegg.jp/list/orthology', KEGGREST:::.matrixParser, ncol = 2)
KEGG.ko_list <- as.data.frame(KEGG.ko_list)
colnames(KEGG.ko_list) <- c('TREM', 'NAME')

KO <- unique(emapper.TERM2GENE$TERM)
KO <- gsub('ko:', '', KO)
emapper.TERM2NAME <- dplyr::filter(KEGG.ko_list, TREM %in% KO)
non_name_term <- setdiff(KO, emapper.TERM2NAME$TREM)
non_name_term <- paste0('ko:', non_name_term)
emapper.TERM2NAME$TREM <- paste0('ko:', emapper.TERM2NAME$TREM)
emapper.TERM2GENE <-  filter(emapper.TERM2GENE, !TERM %in% non_name_term)

db.KO.emapper <- list(TERM2GENE = emapper.TERM2GENE,
                      TERM2NAME = emapper.TERM2NAME)

rm(KO, non_name_term, emapper.TERM2GENE, emapper.TERM2NAME, emapper.Protein2KO)



## build KEGG pathway annotations BY emapper predicted
## --> db.KEGGpathway.emapper

emapper.Protein2Pathway <- emapper_annot %>% 
  dplyr::select(query, KEGG_Pathway) %>% 
  tidyr::drop_na(KEGG_Pathway)

emapper.TERM2GENE <- NULL
for (row in 1:nrow(emapper.Protein2Pathway)) {
  TERM <- str_split(emapper.Protein2Pathway[row, "KEGG_Pathway"], ",", simplify = FALSE)[[1]]  
  GENE <- emapper.Protein2Pathway[row, "query"][[1]]
  tmp <- data.frame(TERM = TERM, GENE = GENE)
  emapper.TERM2GENE <- rbind(emapper.TERM2GENE, tmp)
}; rm(row, GENE, tmp, TERM)
emapper.TERM2GENE <- unique(emapper.TERM2GENE)
emapper.TERM2GENE <- emapper.TERM2GENE[!grepl('^ko', emapper.TERM2GENE$TERM), ]

KEGG.pathway_list <- KEGGREST:::.getUrl('https://rest.kegg.jp/list/pathway', KEGGREST:::.matrixParser, ncol = 2)
KEGG.pathway_list <- as.data.frame(KEGG.pathway_list)
colnames(KEGG.pathway_list) <- c('TREM', 'NAME')

PATHWAY <- unique(emapper.TERM2GENE$TERM)
emapper.TERM2NAME <- dplyr::filter(KEGG.pathway_list, TREM %in% PATHWAY)
non_name_term <- setdiff(PATHWAY, emapper.TERM2NAME$TREM)
emapper.TERM2GENE <-  filter(emapper.TERM2GENE, !TERM %in% non_name_term)

db.KEGGpathway.emapper <- list(TERM2GENE = emapper.TERM2GENE,
                               TERM2NAME = emapper.TERM2NAME)
rm(emapper.TERM2GENE, emapper.TERM2NAME, emapper.Protein2Pathway, PATHWAY, non_name_term)



## build GO annotations BY emapper predicted
## --> db.GO.emapper

emapper.Protein2GO <- emapper_annot %>% 
  dplyr::select(query, GOs) %>% 
  tidyr::drop_na(GOs)

emapper.TERM2GENE <- NULL
for (row in 1:nrow(emapper.Protein2GO)) {
  TERM <- str_split(emapper.Protein2GO[row, "GOs"], ",", simplify = FALSE)[[1]]  
  GENE <- emapper.Protein2GO[row, "query"][[1]]
  tmp <- data.frame(TERM = TERM, GENE = GENE)
  emapper.TERM2GENE <- rbind(emapper.TERM2GENE, tmp)
}; rm(row, GENE, tmp, TERM)
emapper.TERM2GENE <- unique(emapper.TERM2GENE)
emapper.TERM2NAME <- go2term(unique(emapper.TERM2GENE$TERM))

db.GO.emapper <- list(TERM2GENE = emapper.TERM2GENE,
                      TERM2NAME = emapper.TERM2NAME)

rm(emapper.TERM2GENE, emapper.TERM2NAME, emapper.Protein2GO)



## build COG annotations BY emapper predicted
## --> db.COG.emapper

## 项目主页: https://www.ncbi.nlm.nih.gov/research/cog-project/

# emapper_annot 是注释结果数据框
# 包含列：GENE 和 COG_category，其中 COG_category 是以逗号分隔的 COG 分类代码
emapper.Protein2COG <- emapper_annot %>% 
  dplyr::select(query, COG_category) %>% 
  tidyr::drop_na(COG_category)  # 去除 COG_category 为空的行

# 初始化 TERM2GENE
emapper.TERM2GENE <- NULL

# 遍历每一行，将 COG_category 拆分为单独的 COG 分类代码，并与 GENE 映射
for (row in 1:nrow(emapper.Protein2COG)) {
  TERM <- strsplit(emapper.Protein2COG[row, "COG_category"], split = "")[[1]] # 分割多个 COG 分类
  GENE <- emapper.Protein2COG[row, "query"][[1]]  # 提取对应的基因
  tmp <- data.frame(TERM = TERM, GENE = GENE)  # 生成映射数据框
  emapper.TERM2GENE <- rbind(emapper.TERM2GENE, tmp)  # 合并到结果
}; rm(row, GENE, tmp, TERM)

# 去重
emapper.TERM2GENE <- unique(emapper.TERM2GENE)

# 创建 TERM2NAME：COG 分类与分类描述的映射
# TERM2NAME 数据框假设来 https://www.ncbi.nlm.nih.gov/research/cog
COG_descriptions <- data.frame(
  TERM = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S"),
  NAME = c(
    "Translation, ribosomal structure and biogenesis",
    "RNA processing and modification",
    "Transcription",
    "Replication, recombination and repair",
    "Chromatin structure and dynamics",
    "Cell cycle control, cell division, chromosome partitioning",
    "Nuclear structure",
    "Defense mechanisms",
    "Signal transduction mechanisms",
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Cytoskeleton",
    "Extracellular structures",
    "Intracellular trafficking, secretion, and vesicular transport",
    "Posttranslational modification, protein turnover, chaperones",
    "Mobilome: prophages, transposons", 
    "Energy production and conversion",
    "Carbohydrate transport and metabolism",
    "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport and catabolism",
    "General function prediction only",
    "Function unknown"
  )
)

# 过滤 TERM2NAME 中只保留在 TERM2GENE 中的 TERM
emapper.TERM2NAME <- COG_descriptions %>%
  filter(TERM %in% unique(emapper.TERM2GENE$TERM))

# 创建最终的 COG_category 数据库
db.COG.emapper <- list(
  TERM2GENE = emapper.TERM2GENE,
  TERM2NAME = emapper.TERM2NAME
)

rm(emapper.TERM2GENE, emapper.TERM2NAME, emapper.Protein2COG)



## build BRITE annotations BY emapper predicted
## --> db.BRITE.emapper

emapper.Protein2BRITE <- emapper_annot %>% 
  dplyr::select(query, BRITE) %>% 
  tidyr::drop_na(BRITE)  # 去除 BRITE 为空的行

# 初始化 TERM2GENE
emapper.TERM2GENE <- NULL

# 遍历每一行，将 BRITE 分类条目拆分并与 GENE 映射
for (row in 1:nrow(emapper.Protein2BRITE)) {
  TERM <- str_split(emapper.Protein2BRITE[row, "BRITE"], ",", simplify = FALSE)[[1]]  # 分割多个 BRITE 分类
  GENE <- emapper.Protein2BRITE[row, "query"][[1]]  # 提取对应的基因
  tmp <- data.frame(TERM = TERM, GENE = GENE)  # 生成映射数据框
  emapper.TERM2GENE <- rbind(emapper.TERM2GENE, tmp)  # 合并到结果
}; rm(row, GENE, tmp, TERM)

# 去重
emapper.TERM2GENE <- unique(emapper.TERM2GENE)

KEGG.brite_list <- KEGGREST:::.getUrl('https://rest.kegg.jp/list/brite', KEGGREST:::.matrixParser, ncol = 2)
KEGG.brite_list <- as.data.frame(KEGG.brite_list)
colnames(KEGG.brite_list) <- c('TERM', 'NAME')

BRITE <- unique(emapper.TERM2GENE$TERM)
emapper.TERM2NAME <- dplyr::filter(KEGG.brite_list, TERM %in% BRITE)
non_name_term <- setdiff(BRITE, emapper.TERM2NAME$TERM)
emapper.TERM2GENE <-  filter(emapper.TERM2GENE, !TERM %in% non_name_term)

# 创建最终的 BRITE 数据库
db.BRITE.emapper <- list(
  TERM2GENE = emapper.TERM2GENE,
  TERM2NAME = emapper.TERM2NAME
)

rm(emapper.TERM2GENE, emapper.TERM2NAME, emapper.Protein2BRITE, BRITE, non_name_term)


enricher_with_select_db <- function(gene, db, ...) {
  require(clusterProfiler)  
  enricher(gene = gene,
           TERM2GENE = db$TERM2GENE,
           TERM2NAME = db$TERM2NAME,
           ...)
}

rm(KEGG.brite_list, KEGG.ko_list, KEGG.pathway_list, COG_descriptions, readDelimWithDoubleChar)

# save.image('./r_image/enricher_db_panaroo.RData')
# save.image('./r_image/enricher_db_panaroo_pan.RData')

# enricher_with_select_db(
#   gene = emapper_annot$GENE[100:200], 
#   db = db.KO.emapper, 
#   pvalueCutoff = 0.5,
#   qvalueCutoff = 0.5
# ) %>% dotplot()
# 
# 
# enricher.KEGGpathway.emapper.res <- enricher_with_select_db(
#   gene = colnames(coverage.df)[colSums(coverage.df) < 440], 
#   db = db.KEGGpathway.emapper, 
#   pvalueCutoff = 0.1,
#   qvalueCutoff = 0.05
# )
# 
# dotplot(enricher.KEGGpathway.emapper.res)
