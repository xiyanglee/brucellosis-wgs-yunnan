library(tidyverse)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



## Lineage level 2 colors
c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E", 
  "2.1" = "#D0A341", "2.2" = "#E9CA93", "2.3" = "#FEEDB5", 
  "3.1" = "#c9a7b7", "3.2" = "#EBCACB")

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

write.csv(df, "./r_plot_data/figure_2d~df.csv", row.names = FALSE)
