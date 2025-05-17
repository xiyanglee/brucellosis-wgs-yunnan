library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



## Lineage level 2 colors
c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
  "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
  "3.1" = "#c9a7b7", "3.2" = "#EBCACB")

## Lineage level 2a colors
c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3a" = "#9079ad", "1.3b" = "#745399", 
  "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
  "3.1" = "#c9a7b7", "3.2" = "#EBCACB")

metadata.sample %>%
  group_by(collection_year, county_name, lineage_level_2) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collection_year, lineage_level_2) %>%
  summarise(covered_counties = n_distinct(county_name), .groups = "drop") %>% 
  ggplot(aes(x = collection_year, y = covered_counties, color = lineage_level_2)) +
  geom_line(linewidth = 1) +  # 绘制线
  geom_point(size = 3) +  # 绘制点
  labs(x = "Year",
       y = "The number of detected counties each year",
       color = "Lineage") + 
  scale_color_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E"), na.value = "transparent") +
  theme_classic(base_size = 12)

# ggsave('./plot/overview_plot/l2_line~lineage_proportion_of_sampling_time.pdf', width = 3.5, height = 3)

metadata.sample %>%
  group_by(collection_year, county_name, lineage_level_2a) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collection_year, lineage_level_2a) %>%
  summarise(covered_counties = n_distinct(county_name), .groups = "drop") %>% 
  ggplot(aes(x = collection_year, y = covered_counties, color = lineage_level_2a)) +
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(x = "Year",
       y = "The number of detected counties each year",
       color = "Lineage") + 
  scale_color_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3a" = "#9079ad", "1.3b" = "#745399"), na.value = "transparent") +
  theme_classic(base_size = 12)

# ggsave('./plot/overview_plot/l2a_line~lineage_proportion_of_sampling_time.pdf', width = 3.5, height = 3)

metadata.sample %>%
  group_by(collection_year, county_name, lineage_level_2) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collection_year, lineage_level_2) %>%
  summarise(covered_counties = n_distinct(county_name), .groups = "drop") %>%
  group_by(collection_year) %>%
  mutate(total_counties = sum(covered_counties),
         percent = covered_counties / total_counties) %>%
  ggplot(aes(x = factor(collection_year), y = percent, fill = lineage_level_2)) +
  geom_col(position = "stack", color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Year",
    y = "Proportion of counties covered",
    fill = "Lineage"
  ) +
  theme_classic(base_size = 12) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E"), na.value = "transparent")

# ggsave('./plot/overview_plot/l2_bar~lineage_proportion_of_sampling_time.pdf', width = 3.5, height = 3)

library(ggalluvial)

df <- metadata.sample %>%
  group_by(collection_year, county_name, lineage_level_2) %>%
  summarise(.groups = "drop") %>%
  group_by(collection_year, lineage_level_2) %>%
  summarise(covered_counties = n_distinct(county_name), .groups = "drop") %>%
  group_by(collection_year) %>%
  mutate(
    total_counties = sum(covered_counties),
    percent = covered_counties / total_counties
  ) %>%
  ungroup() %>%
  mutate(
    year    = factor(collection_year),
    lineage = lineage_level_2
  )

ggplot(df,
       aes(x = year,
           stratum   = lineage,
           alluvium  = lineage,
           y         = percent,
           fill      = lineage)) +
  # 2.1 流线
  geom_alluvium(
    width   = 0.6,
    alpha   = 0.6,
    knot.pos= 0.5
  ) +
  # 2.2 柱体
  geom_stratum(
    width = 0.6,
    color = "black",
    size  = 0.5
  ) +
  # 2.3 柱体上文字
  geom_text(
    stat  = "stratum",
    aes(label = scales::percent(percent, accuracy = 1)),
    size  = 3, color = "white"
  ) +
  # 2.4 手动配色
  scale_fill_manual(
    values = c(
      "1.1" = "#7EBB8F",
      "1.2" = "#9DD0CF",
      "1.3" = "#886A9E"
    ),
    name = "Lineage"
  ) +
  # 2.5 其余装饰
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Year",
    y = "Proportion of counties covered",
    fill = "Lineage"
  ) + 
  theme_classic(base_size = 12) + 
  coord_cartesian(expand = FALSE) + 
  theme(
    axis.text.x     = element_text(size = 11),
    axis.title      = element_text(size = 13),
    legend.position = "right"
  )

# ggsave('./plot/overview_plot/l2_bar~lineage_proportion_of_sampling_time_type2.pdf', width = 4, height = 3)

metadata.sample %>%
  group_by(collection_year, county_name, lineage_level_2a) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(collection_year, lineage_level_2a) %>%
  summarise(covered_counties = n_distinct(county_name), .groups = "drop") %>%
  group_by(collection_year) %>%
  mutate(total_counties = sum(covered_counties),
         percent = covered_counties / total_counties) %>%
  ggplot(aes(x = factor(collection_year), y = percent, fill = lineage_level_2a)) +
  geom_col(position = "stack", color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Year",
    y = "Proportion of counties covered",
    fill = "Lineage"
  ) +
  theme_classic(base_size = 12) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3a" = "#9079ad", "1.3b" = "#745399"), na.value = "transparent")

# ggsave('./plot/overview_plot/l2a_bar~lineage_proportion_of_sampling_time.pdf', width = 3.5, height = 3)

summarise.patient_type.level2 <- metadata.sample %>% 
  # mutate(patient_type = dplyr::recode(
  #   patient_type,
  #   '其他省' = 'out-of-city cases',
  #   '本省其它地市' = 'out-of-city cases',
  #   '本县区' = 'in-city cases',
  #   '本市其它县区' = 'in-city cases'
  # )) %>%
  mutate(
    patient_type = dplyr::recode(
      patient_type,
      '其他省' = 'out-of-city',
      '本省其它地市' = 'out-of-city',
      '本县区' = 'in-city',
      '本市其它县区' = 'in-city'
    ),
    patient_type = factor(patient_type, levels = c('in-city', 'out-of-city')) # levels = c('in-county', 'out-of-county', 'out-of-city')
  ) %>%
  group_by(patient_type, lineage_level_2) %>% 
  summarise(count = n(), .groups = "drop")
  
ggplot(summarise.patient_type.level2, aes(x = patient_type, y = count, fill = lineage_level_2)) +
  # geom_col(position = "stack", color = "black") + 
  geom_col(position = "fill", color = "black") + 
  labs(x = "Geographic origin of patients", y = "Proportion (%)", fill = "Lineage") +
  theme_classic(base_size = 14) + 
  coord_cartesian(expand = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E"), na.value = "transparent")

# ggsave('./plot/overview_plot/l2_bar~geographic_origin_of_patients.pdf', width = 2.5, height = 3)

summarise.patient_type.level2a <- metadata.sample %>% 
  mutate(
    patient_type = dplyr::recode(
      patient_type,
      '其他省' = 'out-of-city',
      '本省其它地市' = 'out-of-city',
      '本县区' = 'in-city',
      '本市其它县区' = 'in-city'
    ),
    patient_type = factor(patient_type, levels = c('in-city', 'out-of-city')) 
  ) %>%
  group_by(patient_type, lineage_level_2a) %>% 
  summarise(count = n(), .groups = "drop")

ggplot(summarise.patient_type.level2a, aes(x = patient_type, y = count, fill = lineage_level_2a)) +
  # geom_col(position = "stack", color = "black") + 
  geom_col(position = "fill", color = "black") + 
  labs(x = "Geographic origin of patients", y = "Proportion (%)", fill = "Lineage") +
  theme_classic(base_size = 14) + 
  coord_cartesian(expand = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3a" = "#9079ad", "1.3b" = "#745399"), na.value = "transparent")

# ggsave('./plot/overview_plot/l2a_bar~geographic_origin_of_patients.pdf', width = 2.5, height = 3)

metadata.sample %>%
  mutate(profession = case_when(
    profession == "农民"         ~ "Farmer",
    profession == "牧民"         ~ "Herder",
    profession == "散居儿童"     ~ "Non-institutional Child",
    profession == "学生"         ~ "Student",
    profession == "家务及待业"   ~ "Homemaker / Unemployed",
    profession == "离退人员"     ~ "Retired personnel",
    profession == "干部职员"     ~ "Office staff",
    profession == "其他"         ~ "Unknown / Other",
    profession == "不详"         ~ "Unknown / Other",
    TRUE                         ~ NA_character_
  )) %>% 
  group_by(profession, lineage_level_2) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(aes(x = lineage_level_2, y = count, fill = profession)) +
  geom_col(position = "stack", color = "black") +
  # geom_col(position = "fill", color = "black") + 
  labs(x = "Lineage", y = "Case count", fill = "Profession") +
  scale_fill_manual(
    values = c(
      "Farmer" =                 "#82ae46",
      "Herder" =                 "#757cbb",
      "Child" =                  "#d77f66",
      "Student" =                "#da9233",
      "Homemaker / Unemployed" = "#6f94cd",
      "Retired personnel" =      "#5aa4ae",
      "Office staff" =           "#ec6d51",
      "Unknown / Other" =        "#e7e7eb"
    )
  ) + 
  theme_minimal(base_size = 14)

# ggsave('./plot/overview_plot/l2_bar~case_count_of_profession.pdf', width = 5, height = 3)

metadata.sample %>%
  mutate(profession = case_when(
    profession == "农民"         ~ "Farmer",
    profession == "牧民"         ~ "Herder",
    profession == "散居儿童"     ~ "Non-institutional Child",
    profession == "学生"         ~ "Student",
    profession == "家务及待业"   ~ "Homemaker / Unemployed",
    profession == "离退人员"     ~ "Retired personnel",
    profession == "干部职员"     ~ "Office staff",
    profession == "其他"         ~ "Unknown / Other",
    profession == "不详"         ~ "Unknown / Other",
    TRUE                         ~ NA_character_
  )) %>% 
  group_by(profession, lineage_level_2a) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(aes(x = lineage_level_2a, y = count, fill = profession)) +
  geom_col(position = "stack", color = "black") +
  # geom_col(position = "fill", color = "black") + 
  labs(x = "Lineage", y = "Case count", fill = "Profession") +
  scale_fill_manual(
    values = c(
      "Farmer" =                 "#82ae46",
      "Herder" =                 "#757cbb",
      "Child" =                  "#d77f66",
      "Student" =                "#da9233",
      "Homemaker / Unemployed" = "#6f94cd",
      "Retired personnel" =      "#5aa4ae",
      "Office staff" =           "#ec6d51",
      "Unknown / Other" =        "#e7e7eb"
    )
  ) + 
  theme_minimal(base_size = 14)

# ggsave('./plot/overview_plot/l2a_bar~case_count_of_profession.pdf', width = 5.5, height = 3)
