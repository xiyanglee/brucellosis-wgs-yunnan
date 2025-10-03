library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')



pangenes.roary <- read.csv('./data/roary_result/gene_presence_absence.csv', header = TRUE, stringsAsFactors = FALSE)
allgenes.roary <- pangenes.roary %>%
  dplyr::rename(roary_gene = Gene) %>%  # 重命名 Gene 列
  pivot_longer(
    # cols = matches("^(X|GCA_)"),  # 筛选以 X 或 GCA_ 开头的列
    cols = matches("^X"), 
    names_to = "sample_name",     # 新列名：原列名变为 sample_name
    values_to = "prokka_gene"     # 原值变为 prokka_gene
  ) %>% 
  dplyr::select(sample_name, prokka_gene, roary_gene) %>% 
  separate_rows(prokka_gene, sep = "\t")
allgenes.roary <- allgenes.roary[allgenes.roary$prokka_gene != "", ]

pangenes.panaroo <- read.csv('./data/panaroo_result/gene_presence_absence_roary.csv', header = TRUE, stringsAsFactors = FALSE)
allgenes.panaroo <- pangenes.panaroo %>%
  dplyr::rename(panaroo_gene = Gene) %>%  # 重命名 Gene 列
  pivot_longer(
    cols = matches("^X"), 
    names_to = "sample_name", 
    values_to = "prokka_gene" 
  ) %>% 
  dplyr::select(sample_name, prokka_gene, panaroo_gene) %>% 
  # 拆解 prokka_gene 字段按分号分行
  separate_rows(prokka_gene, sep = ";") %>% 
  dplyr::filter(prokka_gene != "")

rm(pangenes.roary, pangenes.panaroo)
gc()


length(setdiff(allgenes.panaroo$prokka_gene, allgenes.roary$prokka_gene))
# [1] 2928
head(setdiff(allgenes.panaroo$prokka_gene, allgenes.roary$prokka_gene), 10)
# [1] "2_refound_47"    "5_refound_142"   "2_refound_46"    "5_refound_141"   "1_refound_11"    "14_refound_421"  "23_refound_691"  "29_refound_823" 
# [9] "72_refound_1967" "97_refound_2718"

length(setdiff(allgenes.roary$prokka_gene, allgenes.panaroo$prokka_gene))
# [1] 3248
head(setdiff(allgenes.roary$prokka_gene, allgenes.panaroo$prokka_gene), 10)
# [1] "PROKKATAG-201915_02540" "PROKKATAG-201917_02904" "PROKKATAG-201924_02807" "PROKKATAG-201928_02533" "PROKKATAG-201931_02809" "PROKKATAG-201937_02537"
# [7] "PROKKATAG-202003_02141" "PROKKATAG-202013_02544" "PROKKATAG-202018_02892" "PROKKATAG-202019_02535"

length(intersect(allgenes.panaroo$prokka_gene, allgenes.roary$prokka_gene))
# [1] 318958


allgenes <- full_join(allgenes.panaroo, allgenes.roary, by = c("sample_name", "prokka_gene"))
rm(allgenes.panaroo, allgenes.roary)
gc()

allgenes$sample_name <- gsub("^X", "", allgenes$sample_name)

saveRDS(allgenes, './r_image/panaroo_roary_bridge.RData')

