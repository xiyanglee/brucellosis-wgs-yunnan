load('./r_image/metadata.RData')
load('./r_image/cluster_tree.RData')

metadata_tree <- left_join(group_data_2, metadata, by = c('tip' = 'assembly_accession'))

report_card.yunnan <- readxl::read_xlsx('../../metadata/云南布病病例数据2006-2022-analyses.xlsx', sheet = '2006-2022年病例', col_types = "text")
report_card.yunnan <- report_card.yunnan %>% 
  dplyr::select(现住详细地址, 发病日期) %>%
  mutate(发病日期 = as.Date(as.integer(发病日期), origin = "1899-12-30"), 
         city = substr(现住详细地址, 4, 5), 
         city = recode(
           city, '西双' = '西双版纳'
         ), 
         city_full = recode(
           city, 
           '丽江' = '丽江市', '临沧' = '临沧市', '保山' = '保山市',  '大理' = '大理白族自治州',  '德宏' = '德宏傣族景颇族自治州',  '怒江' = '怒江傈僳族自治州',  
           '文山' = '文山壮族苗族自治州',  '昆明' = '昆明市',  '昭通' = '昭通市',  '普洱' = '普洱市',  '曲靖' = '曲靖市',  
           '楚雄' = '楚雄彝族自治州',  '玉溪' = '玉溪市', '红河' = '红河哈尼族彝族自治州', '西双版纳' = '西双版纳傣族自治州', '迪庆' = '迪庆藏族自治州', 
         ), 
         city_en = recode(
           city, 
           '丽江' = 'LiJiang', '临沧' = 'LinCang', '保山' = 'BaoShan',  '大理' = 'DaLi',  '德宏' = 'DeHong',  '怒江' = 'NuJiang',  
           '文山' = 'WenShan',  '昆明' = 'KunMing',  '昭通' = 'ZhaoTong',  '普洱' = 'PuEr',  '曲靖' = 'QuJing',  
           '楚雄' = 'ChuXiong',  '玉溪' = 'YuXi', '红河' = 'HongHe', '西双版纳' = 'XiShuangBanNa', '迪庆' = 'DiQing', 
           .default = 'Others'
         )) %>% 
  dplyr::filter(city_en != 'Others') %>% 
  select(-c(现住详细地址, city))
colnames(report_card.yunnan) <- c('onset_date', 'city_full', 'city_en')
report_card.yunnan <- report_card.yunnan %>% 
  mutate(
    collection_year  = year(onset_date), 
    collection_month = month(onset_date)
  ) %>% 
  select(-onset_date) %>% 
  group_by(collection_year, collection_month) %>% 
  summarize(incidence = n(), .groups = 'drop') %>% 
  ungroup()

report_card.yunnan <- report_card.yunnan %>% 
  mutate(collection_date = as.numeric(collection_year) + collection_month / 12) %>% 
  dplyr::filter(collection_year >= 2019)

summary_tree <- metadata_tree %>% 
  filter(geographic_country == 'Yunnan') %>% 
  mutate(collection_date = as.numeric(collection_year) + collection_month / 12) %>% 
  left_join(., dplyr::select(report_card.yunnan, collection_date, incidence), by = 'collection_date')

# report_card.yunnan <- report_card.yunnan %>% 
#   mutate(collection_date = as.Date(paste(collection_year, collection_month, '1', sep = "-"), format = "%Y-%m-%d")) %>% 
#   dplyr::filter(collection_year >= 2019)
# 
# summary_tree <- metadata_tree %>% 
#   filter(geographic_country == 'Yunnan') %>% 
#   group_by(collection_year, collection_month) %>% 
#   summarize(incidence = n(), .groups = 'drop') %>% 
#   mutate(collection_date = as.Date(paste(collection_year, collection_month, '1', sep = "-"), format = "%Y-%m-%d"))

ggplot() + 
  geom_line(data = report_card.yunnan, aes(x = collection_date, y = incidence)) + 
  # geom_line(data = summary_tree, aes(x = collection_date, y = incidence)) + 
  geom_segment(data = summary_tree, aes(x = collection_date, xend = collection_date, y = 0, yend = incidence), linewidth = 1, alpha = 0.2) + 
  theme_classic()



# sampling_patients.yunnan <- readxl::read_xlsx('../../metadata/云南布病病例数据2006-2022-analyses.xlsx', sheet = '菌株测序', col_types = "text") %>% 
#   dplyr::select(发病日期) %>% 
#   mutate(发病日期 = as.Date(as.integer(发病日期), origin = "1899-12-30")) %>% 
#   rename(onset_date = `发病日期`)

sampling_patients.yunnan <- dplyr::filter(metadata_tree, geographic_country == 'Yunnan') %>% 
  mutate(collection_date = as.Date(paste(collection_year, collection_month, collection_day, sep = "-"), format = "%Y-%m-%d"))


report_card.yunnan <- report_card.yunnan %>%
  mutate(collection_date = as.Date(paste(collection_year, collection_month, '1', sep = "-"), format = "%Y-%m-%d")) %>%
  dplyr::filter(collection_year >= 2019)

ggplot() + 
  geom_segment(data = sampling_patients.yunnan, aes(x = collection_date, xend = collection_date, y = 0, yend = 140, color = group_2), linewidth = 0.5, alpha = 1) + 
  geom_line(data = report_card.yunnan, aes(x = collection_date, y = incidence), linewidth = 1) + 
  geom_ribbon(data = report_card.yunnan, aes(x = collection_date, ymin = pmin(incidence, 140), ymax = 140), fill = "white", alpha = 1) + 
  # geom_line(data = summary_tree, aes(x = collection_date, y = incidence)) + 
  theme_classic(base_size = 14) + 
  scale_color_manual(values = c("1.1" = "#FBCBD2", "1.2" = "#ee7959", "1.3" = "#ab1d22", 
                                "2.1" = "#aed0ee", "2.2" = "#2ca9e1", "2.3" = "#2e59a7", 
                                "3.1" = "#c0d695", "3.2" = "#4f6f46"), 
                     na.value = "transparent") + 
  labs(
    x = "Date",
    y = "Incidence",
    color = "Lineage",
    title = NULL
  )

# ggsave('./plot/overview_plot/sampling_time.pdf', width = 4.5, height = 3)


table(sampling_patients.yunnan$group_4, sampling_patients.yunnan$collection_year)

sampling_patients.yunnan %>% 
  dplyr::select(group_2, collection_year) %>% 
  na.omit() %>% 
  group_by(group_2, collection_year) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(collection_year) %>%
  mutate(proportion = count / sum(count)) %>% 
  ggplot(aes(x = collection_year, y = proportion, fill = group_2)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Year",
    y = "Proportion",
    fill = "Lineage",
    title = NULL
  ) +
  theme_classic(base_size = 10) + 
  scale_fill_manual(values = c("1.1" = "#FBCBD2", "1.2" = "#ee7959", "1.3" = "#ab1d22", 
                               "2.1" = "#aed0ee", "2.2" = "#2ca9e1", "2.3" = "#2e59a7", 
                               "3.1" = "#c0d695", "3.2" = "#4f6f46"), 
                    na.value = "transparent")

# ggsave('./plot/overview_plot/lineage_proportion_of_sampling_time.pdf', width = 3.5, height = 3)

