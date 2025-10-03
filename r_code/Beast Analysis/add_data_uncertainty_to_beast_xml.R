library(xml2)

setwd('/Users/xiyangli/Lab/Project/Brucella_WGS_Yunnan_2024/')



metadata <- read.csv('./pipelines/denovo_ssembly/data/beast_gubbins_contig/metadata.csv')

# metadata_loc <- dplyr::select(metadata, c(assembly_accession, area))
# colnames(metadata_loc) <- c('traits', 'state')
# metadata_loc$state <- gsub(' ', '_', metadata_loc$state)
# write.table(metadata_loc, './pipelines/denovo_ssembly/data/beast_gubbins_contig/location.txt', quote = FALSE, sep = '\t', row.names = FALSE)

ggplot(metadata, aes(x = area, y = collection_date)) + 
  geom_point() + 
  coord_cartesian(ylim = c(2015, 2023))

metadata <- metadata %>% 
  mutate(trait_area = case_when(
    area == 'East Asia' ~ 'East Asia', 
    area %in% c('Southeast Asia', 'South Asia') ~ 'Southeast Asia', 
    area == 'North Asia' ~ 'North Asia', 
    area %in% c('Western Asia', 'Central Asia') ~ 'Western and Central Asia', 
    .default = continent
  )) 

ggplot(metadata, aes(x = trait_area, y = collection_date)) + 
  geom_point() + 
  coord_cartesian(ylim = c(2015, 2024))

metadata_loc <- dplyr::select(metadata, c(assembly_accession, trait_area))
colnames(metadata_loc) <- c('traits', 'state')
metadata_loc$state <- gsub(' ', '_', metadata_loc$state)
write.table(metadata_loc, './pipelines/denovo_ssembly/data/beast_gubbins_contig/location.txt', quote = FALSE, sep = '\t', row.names = FALSE)



# 读取 XML 文件
# xml_file <- read_xml('./pipelines/denovo_ssembly/data/beast/core_gene_snp_filtered.xml')
# metadata <- read.csv('./pipelines/denovo_ssembly/data/beast/metadata.csv')
xml_file <- read_xml('./pipelines/denovo_ssembly/data/beast_gubbins_contig/beast_GTRIG_clean_core_filtered.xml')

# 遍历所有 taxon 节点
taxon_nodes <- xml_find_all(xml_file, ".//taxon")

for (taxon in taxon_nodes) {
  # 获取 taxon 的 id 属性
  taxon_id <- xml_attr(taxon, "id")
  
  if (!is.na(taxon_id)) {
    date_node <- xml_find_first(taxon, "./date")
    
    row_id <- match(taxon_id, metadata$assembly_accession)
    uncertainty_value <- metadata[row_id, 'uncertainty']
    date_value <- metadata[row_id, 'collection_date']
    
    xml_set_attr(date_node, "value", date_value)
    xml_set_attr(date_node, "uncertainty", uncertainty_value)
  }
}

write_xml(xml_file, "./pipelines/denovo_ssembly/data/beast_gubbins_contig/beast_GTRIG_clean_core_filtered_with_time.xml")
