# install.packages("rentrez", lib = "/public/r_share_library/4.1")
# install.packages("biomartr", lib = "/public/r_share_library/4.1")
library(tidyverse)
library(xml2)
library(purrr)
library(rentrez)
# library(biomartr)
library(progress)

load('./r_image/metadata_v3.RData')



# 设定 assembly_accession 列表
accessions <- metadata.raw.genebank$assembly_accession

# 设置进度条
pb <- progress_bar$new(
  format = "  downloading [:bar] :current/:total (:percent) eta: :eta",
  total = length(accessions), clear = FALSE, width = 60
)

# 单个 accession 的下载函数（最多重试 3 次）
fetch_one_metadata_retry <- function(acc, max_attempts = 3) {
  attempt <- 1
  while (attempt <= max_attempts) {
    result <- tryCatch({
      uid <- entrez_search(db = "assembly", term = acc)$ids[1]
      if (is.null(uid) || is.na(uid)) return(NULL)
      
      summary_info <- entrez_summary(db = "assembly", id = uid)
      return(as.list(summary_info))
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(result)) return(result)
    Sys.sleep(1)  # 间隔时间防止过载
    attempt <- attempt + 1
  }
  return(NULL)
}

# 主循环：带进度条
results_list <- vector("list", length(accessions))
for (i in seq_along(accessions)) {
  results_list[[i]] <- fetch_one_metadata_retry(accessions[i])
  pb$tick()
}; rm(i, pb)

# 过滤无效结果
metalist.entrez <- plyr::compact(metalist.entrez)

# 转为 data.frame（字段可能不一致时使用 bind_rows）
# metadata_df <- bind_rows(metalist.entrez)



metadata.biosample <- map_dfr(metalist.entrez, function(x) {
  tibble(
    assembly_accession = x$assemblyaccession,
    biosample_id = x$biosampleid
  )
})

# 初始化进度条
pb <- progress_bar$new(
  format = "Retrieving [:bar] :current/:total (:percent) eta: :eta",
  total = nrow(metadata.biosample), clear = FALSE, width = 60
)

# 定义带重试的获取函数
get_biosample_metadata <- function(biosample_id, max_attempts = 3) {
  attempt <- 1
  while (attempt <= max_attempts) {
    result <- tryCatch({
      summary_info <- entrez_summary(db = "biosample", id = biosample_id)
      as.list(summary_info)
    }, error = function(e) NULL)
    
    if (!is.null(result)) return(result)
    Sys.sleep(1)
    attempt <- attempt + 1
  }
  return(NULL)
}

metalist.biosample <- vector("list", length = nrow(metadata.biosample))
for (i in seq_len(nrow(metadata.biosample))) {
  metalist.biosample[[i]] <- get_biosample_metadata(metadata.biosample$biosample_id[i])
  pb$tick()
}; rm(i, pb)


# 提取函数
extract_biosample_info <- function(x) {
  # 解析 XML
  xml <- read_xml(x$sampledata)
  
  # 提取两个属性（可能不存在，需要做 NA 处理）
  collection_date <- xml %>%
    xml_find_first(".//Attribute[@attribute_name='collection_date']") %>%
    xml_text()
  
  geo_loc_name <- xml %>%
    xml_find_first(".//Attribute[@attribute_name='geo_loc_name']") %>%
    xml_text()
  
  # 返回一行数据
  tibble(
    uid = x$uid,
    date = x$date,
    collection_date = ifelse(length(collection_date) == 0, NA, collection_date),
    geographic_location = ifelse(length(geo_loc_name) == 0, NA, geo_loc_name)
  )
}

# 应用到所有 biosample
metadata.biosample.info <- map_dfr(metalist.biosample, extract_biosample_info)
metadata.biosample.info <- metadata.biosample.info %>%
  mutate(across(everything(), ~ifelse(. %in% c("Missing", "not applicable", "unknown"), NA, .)))

sum(!is.na(metadata.biosample.info$collection_date))
# [1] 423
sum(!is.na(metadata.raw.genebank$collection_date))
# [1] 608

sum(!is.na(metadata.biosample.info$geographic_location))
# [1] 425
sum(!is.na(metadata.raw.genebank$geographic_location))
# [1] 711

# save.image('./r_image/metadata_v3_with_entrez.RData')
