library(tidyverse)
library(sf)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')



# 由于第一次读取表格没有加 quote = NULL，导致行数不全
# metadata <- read.table('./source/assembly_genbank/ncbi_assembly_brucella_melitensis_metadata.tsv', sep = '\t', quote = NULL, header = TRUE)
# metadata <- dplyr::select(metadata, `Assembly.Accession`, `Assembly.BioSample.Accession`, collection_date, geographic_location)
# colnames(metadata) <- c('assembly_accession', 'biosample_accession', 'collection_date', 'geographic_location')
# 
# 第一次导出的表填充后合并到全表中，再次导出
# metadata_imputed <- read.table('./source/assembly_genbank/ncbi_assembly_brucella_melitensis_metadata_imputed.tsv', sep = '\t', quote = "\"", header = TRUE)
# metadata <- metadata %>%
#   left_join(metadata_imputed, by = c("assembly_accession", "biosample_accession"), suffix = c("_b", "_a")) %>%      # 合并表格，保留后缀
#   mutate(geographic_location = coalesce(geographic_location_a, geographic_location_b), 
#          collection_date     = coalesce(collection_date_a, collection_date_b)) %>%        # 用 a 的值更新 b 的值
#   select(assembly_accession, biosample_accession, collection_date, geographic_location) 
# 
# write.table(metadata, './source/assembly_genbank/ncbi_assembly_brucella_melitensis_metadata_filter.tsv', sep = '\t', row.names = FALSE)



load('./r_image/clustering_analysis_fastBAPS.RData')
rm(list = setdiff(ls(), c("raxml_tree", "model", "result")))

metadata <- read.table('./source/assembly_genbank/ncbi_assembly_brucella_melitensis_metadata_imputed.tsv', sep = '\t', quote = "\"", header = TRUE)

metadata <- metadata %>% 
  separate(collection_date, into = c("collection_year", "collection_month", "collection_day"), sep = "-", remove = FALSE) %>% 
  separate(geographic_location, into = c("geographic_country", "geographic_province"), sep = ":\\s*", remove = FALSE) # \\s* 表示匹配零个或多个空格

metadata <- metadata %>% 
  mutate(
    continent = case_when(
      geographic_country %in% c('Afghanistan', 'China', 'Cyprus', 'India', 'Iran', 'Iraq', 'Israel', 'Jordan', 'Kuwait', 'Malaysia', 'Pakistan', 
                                'Russia', 'Saudi Arabia', 'Syria', 'Thailand', 'Turkey', 'Turkmenistan') ~ 'Asia', 
      geographic_country %in% c('Albania', 'Austria', 'Belgium', 'Bulgaria', 'Croatia', 'France', 'Georgia', 'Italy', 'Kosovo', 'Malta', 'Norway', 
                                'Portugal', 'Serbia', 'United Kingdom') ~ 'Europe', 
      geographic_country %in% c('Argentina', 'Canada', 'USA') ~ 'America', 
      geographic_country %in% c('Egypt', 'Morocco', 'Nigeria', 'Somalia', 'South Africa', 'Sudan', 'Zimbabwe') ~ 'Africa'
    ), 
    area = case_when(
      geographic_country %in% c('Afghanistan', 'Turkmenistan') ~ 'Central Asia', 
      geographic_country %in% c('China') ~ 'East Asia', 
      geographic_country %in% c('Malaysia', 'Thailand') ~ 'Southeast Asia', 
      geographic_country %in% c('Cyprus', 'Iran', 'Iraq', 'Israel', 'Jordan', 'Kuwait', 'Saudi Arabia', 'Syria', 'Turkey') ~ 'Western Asia', 
      geographic_country %in% c('Russia') ~ 'North Asia', 
      geographic_country %in% c('India', 'Pakistan') ~ 'South Asia', 
      geographic_country %in% c('Belgium', 'France', 'Portugal', 'United Kingdom') ~ 'Western Europe', 
      geographic_country %in% c('Georgia') ~ 'Eastern Europe', 
      geographic_country %in% c('Albania', 'Bulgaria', 'Croatia', 'Italy', 'Kosovo', 'Malta', 'Serbia') ~ 'Southern Europe', 
      geographic_country %in% c('Austria') ~ 'Central Europe', 
      geographic_country %in% c('Norway') ~ 'North Europe', 
      geographic_country %in% c('Argentina') ~ 'South America', 
      geographic_country %in% c('Canada', 'USA') ~ 'North America', 
      geographic_country %in% c('Egypt', 'Morocco') ~ 'North Africa', 
      geographic_country %in% c('Nigeria') ~ 'Western Africa', 
      geographic_country %in% c('Somalia', 'Sudan') ~ 'Eastern Africa', 
      geographic_country %in% c('South Africa', 'Zimbabwe') ~ 'Southern Africa'
    )
  )

# 将 metadata 中所有值为 "" 的单元格替换为 NA
metadata <- metadata %>%
  mutate_all(~ na_if(., ""))
#  将 年月日 转换为数值格式
metadata <- metadata %>%
  mutate(across(c(collection_year, collection_month, collection_day), as.integer))

metadata.raw.genebank <- metadata
rm(metadata)



metadata.raw.sample <- readxl::read_xlsx('../../metadata/研究菌株.xlsx', sheet = 'Sheet2')
metadata.raw.sample$`菌株编号` <- as.character(metadata.raw.sample$`菌株编号`)

metadata.sample <- metadata.raw.sample %>% 
  dplyr::select(`菌株编号`, `发病日期`) %>% 
  dplyr::rename(
    assembly_accession = `菌株编号`,
    collection_date = `发病日期`
  ) %>% 
  dplyr::select(assembly_accession, collection_date) %>% 
  mutate(collection_year = year(collection_date), 
         collection_month = month(collection_date), 
         collection_day = day(collection_date)) %>% 
  mutate(geographic_country = 'China', 
         continent = 'Asia', 
         area = 'East Asia')

metadata.sample$collection_date <- as.character(metadata.sample$collection_date)
metadata.sample$this_study <- 'Yunnan'

metadata <- dplyr::bind_rows(metadata.raw.genebank, metadata.sample)
metadata <- metadata[, colnames(metadata.sample)]


partition.merge <- readRDS('./r_image/clustering_analysis_designating_lineages.RDS')

metadata <- metadata %>% 
  left_join(partition.merge, by = c("assembly_accession" = "tip"))

metadata.phylo <- metadata[!is.na(metadata$lineage_level_1), ]
metadata.sample <- metadata.sample %>% 
  left_join(partition.merge, by = c("assembly_accession" = "tip"))



## 提取所有病例信息
metadata.cases <- readxl::read_xlsx('../..//metadata/云南布病病例数据2006-2022-analyses.xlsx', sheet = '2006-2022年病例')
# colnames(metadata.cases)
metadata.cases <- metadata.cases %>% 
  dplyr::select(`卡片编号`, `有效证件号`, `性别`, `出生日期`, `年龄`, `人群分类`, `病人属于`, `现住详细地址`, `现住地址国标`, `发病日期`, `报告单位地区编码`, `审核状态`) %>% 
  dplyr::filter(`审核状态` != '已删除卡') %>% 
  dplyr::select(-`审核状态`) %>% 
  dplyr::rename(
    ID = `卡片编号`, 
    card_ID = `有效证件号`, 
    gender = `性别`, 
    birth_date = `出生日期`, 
    age = `年龄`, 
    address = `现住详细地址`, 
    address_code = `现住地址国标`, 
    onset_date = `发病日期`, 
    patient_type = `病人属于`, 
    profession = `人群分类`, 
    reporting_agency_region_code = `报告单位地区编码`
  )

str(metadata.cases)
metadata.cases$age <- as.integer(gsub("岁", "", metadata.cases$age))
# 年龄为 6月、4月
metadata.cases$age[is.na(metadata.cases$age)] <- 0

## 提取测序的病例信息
metadata.cases.sequencing <- readxl::read_xlsx('../..//metadata/研究菌株.xlsx', sheet = 'Sheet2')
metadata.cases.sequencing <- metadata.cases.sequencing %>% 
  dplyr::select(`有效证件号`, `菌株编号`, `发病日期`) %>% 
  dplyr::rename(
    card_ID = `有效证件号`, 
    tip = `菌株编号`, 
    onset_date = `发病日期`
  )
metadata.cases.sequencing$tip <- gsub(".*[（(]\\s*(\\d+)\\s*[）)].*", "\\1", metadata.cases.sequencing$tip)
# str(sequencing_data)

metadata.cases <- left_join(metadata.cases, metadata.cases.sequencing, by = c('card_ID', 'onset_date')) %>% 
  left_join(., partition.merge, by = 'tip')

metadata.cases.sequencing <- metadata.cases[!is.na(metadata.cases$tip), ]



# 获取两个数据框中同名的非键列
common_cols <- intersect(
  setdiff(names(metadata.sample), "assembly_accession"),
  setdiff(names(metadata.cases.sequencing), "tip")
)

# 从 metadata.cases.sequencing 中移除这些同名列
metadata.cases.sequencing.filtered <- metadata.cases.sequencing %>%
  select(-all_of(common_cols))

# 执行合并
metadata.sample <- left_join(
  metadata.sample,
  metadata.cases.sequencing.filtered,
  by = c("assembly_accession" = "tip")
)

rm(common_cols, metadata.cases.sequencing.filtered, partition.merge)



# 1. 读取菌株采样数据和 geojson 文件
seq_loc <- read.csv('../../metadata/Yunnan_sequenced_patient_location.csv')
seq_loc$input_address <- NULL
colnames(seq_loc) <- c('tip', 'lng', 'lat')
seq_loc$tip <- as.character(seq_loc$tip)

metadata.sample <- left_join(metadata.sample, seq_loc, by = c('assembly_accession' = 'tip'))

# 读取云南省县级 geojson 文件
sf_county <- read_sf('../../metadata/云南省_县.geojson')

# 2. 将菌株采样点转换为 sf 对象（空间点数据）
seq_sf <- st_as_sf(metadata.sample, coords = c("lng", "lat"), crs = st_crs(sf_county))

# 3. 通过空间叠加将菌株分配到具体的县
# 确保 seq_sf 和 sf_county 坐标系一致
seq_sf    <- st_transform(seq_sf, st_crs(sf_county))
sf_county <- st_transform(sf_county, st_crs(sf_county))

# 空间叠加，分配每个点到县
seq_with_county <- st_join(seq_sf, sf_county, join = st_within)

metadata.sample <- left_join(
  metadata.sample, 
  seq_with_county %>% 
    as.data.frame() %>% 
    dplyr::select(assembly_accession, name, gb) %>% 
    dplyr::rename(county_name = name, 
                  county_gb = gb), 
  by = 'assembly_accession'
)

rm(seq_with_county, seq_sf, sf_county, seq_loc)

# metadata.phylo 中添加省份信息 (包括 Yunnan)
metadata.phylo <- left_join(metadata.phylo, metadata.raw.genebank[, c('assembly_accession', 'geographic_province')], by = 'assembly_accession')
metadata.phylo <- metadata.phylo %>% 
  mutate(
    geographic_province = replace_na(geographic_province, ""),
    this_study = replace_na(this_study, ""),
    geographic_province = paste0(geographic_province, this_study)
  )


README <- function() {
  txt <- 
  "
  分析中使用的信息:
  metadata        为总的样本信息，包括公共序列和云南测序样本信息
  metadata.phylo  为已分型的公共序列样本信息
  metadata.sample 为已分型的云南测序样本信息 (包括 metadata.cases.sequencing 中的信息)
  
  仅作备份的信息:
  metadata.raw.genebank     为公共序列样本的原始信息
  metadata.raw.sample       为云南测序样本的原始信息
  metadata.cases            所有患者的基本信息 (原 cases_data)
  metadata.cases.sequencing 测序患者的基本信息，metadata.cases 的子集 (原 sequencing_data)
  
  metadata.sample 各列含义: 
    assembly_accession = `菌株编号`,
    collection_date = `发病日期`
    ID = `卡片编号`, 
    card_ID = `有效证件号`, 
    gender = `性别`, 
    birth_date = `出生日期`, 
    age = `年龄`, 
    address = `现住详细地址`, 
    address_code = `现住地址国标`, 
    onset_date = `发病日期`, 
    patient_type = `病人属于`, 
    profession = `人群分类`, 
    reporting_agency_region_code = `报告单位地区编码`
  
  注释:
  - 首先读取 GeneBank metadata 形式的: metadata.raw.genebank
  - 云南采样样本的 metadata 也按照 GeneBank metadata 的形式生成: metadata.raw.sample
  - 二者合并为 metadata，并合并 Lineage 注释信息
  - metadata 拆分为 metadata.phylo 和 metadata.sample
  - 读取云南省所有布病患者信息 metadata.cases，合并 Lineage 注释信息
  - 从 metadata.cases 提取测序样本子集 metadata.cases.sequencing
  - metadata.sample 合并 metadata.cases.sequencing 中的新增信息
  - metadata.sample 新增患者住址坐标，并分配到县
  "
  cat(txt)
}


README()

# save.image('./r_image/metadata_v3.RData')
