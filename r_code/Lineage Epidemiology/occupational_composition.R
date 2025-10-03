library(tidyverse)

setwd('/Users/xiyangli/Lab/Project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



df_clean <- metadata.sample %>%
  mutate(profession = case_when(
    profession == "农民"         ~ "a)",
    profession == "牧民"         ~ "a)",
    profession == "散居儿童"     ~ "c)",
    profession == "教师"         ~ "c)",
    profession == "学生"         ~ "c)",
    profession == "医务人员"     ~ "c)",
    profession == "商业服务"     ~ "c)",
    profession == "工人"         ~ "c)",
    profession == "幼托儿童"     ~ "c)",
    profession == "民工"         ~ "c)",
    profession == "家务及待业"   ~ "c)",
    profession == "离退人员"     ~ "c)",
    profession == "干部职员"     ~ "c)",
    profession == "餐饮食品业"   ~ "c)",
    profession %in% c("其他","不详") ~ "b)",
    TRUE                         ~ NA_character_
  )) %>% 
  filter(lineage_level_2 %in% c("1.1", "1.3")) %>%
  filter(!is.na(lineage_level_2), !is.na(profession))

# 列联表：行=lineage_level_2，列=profession
tab <- df_clean %>%
  count(lineage_level_2, profession) %>%
  pivot_wider(names_from = profession, values_from = n, values_fill = 0) %>%
  as.data.frame()

rownames(tab) <- tab$lineage_level_2
tab$lineage_level_2 <- NULL
tab <- as.matrix(tab)

fisher.test(tab)
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  tab
# p-value = 0.04771
# alternative hypothesis: two.sided


