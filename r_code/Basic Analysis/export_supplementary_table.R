# install.packages("openxlsx", lib = "/public/r_share_library/4.1")

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata.RData')
load('./r_image/cluster_tree.RData')

left_join(metadata.genebank, group_data_2, by = c("assembly_accession" = "tip")) %>% 
  select(assembly_accession, collection_date, geographic_country, geographic_province, continent, area, group_1, group_2) %>% 
  rename(`Assembly Accession` = assembly_accession, 
         `Colleation Date` = collection_date, 
         `Colleation Country` = geographic_country, 
         `Colleation Province` = geographic_province, 
         `Continent` = continent, 
         `Area` = area, 
         `Lineage` = group_1, 
         `Sub-lineage` = group_2) %>% 
  openxlsx::write.xlsx("./result/supplementary_table/supplementary_table_1.xlsx", sheetName = "Supplementary Table 1")

left_join(metadata.sample, group_data_2, by = c("assembly_accession" = "tip")) %>% 
  mutate(geographic_province = 'Yunnan', 
         geographic_country = 'China', 
         continent = 'Asia', 
         area = 'East Asia') %>% 
  select(assembly_accession, collection_date, geographic_country, geographic_province, continent, area, group_1, group_2) %>% 
  rename(`Assembly Accession` = assembly_accession, 
         `Colleation Date` = collection_date, 
         `Colleation Country` = geographic_country, 
         `Colleation Province` = geographic_province, 
         `Continent` = continent, 
         `Area` = area, 
         `Lineage` = group_1, 
         `Sub-lineage` = group_2) %>% 
  openxlsx::write.xlsx("./result/supplementary_table/supplementary_table_2.xlsx", sheetName = "Supplementary Table 2")
