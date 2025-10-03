library(tidyverse)
library(ggpmisc)

setwd('/Users/xiyangli/Lab/Project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/distance_dynamic_random_points.RData')

summarize_data <- metadata.sample %>% 
  group_by(collection_year, lineage_level_2) %>% 
  summarize(sample_nums = n(), 
            .groups = "drop") %>% 
  group_by(collection_year) %>% 
  mutate(yearly_total_samples = sum(sample_nums)) %>% 
  ungroup() 

# c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E")

ggplot(data = dplyr::filter(summarize_data, lineage_level_2 == "1.3"), aes(x = collection_year, y = sample_nums, weight = yearly_total_samples)) + 
  geom_smooth(method = "lm", se = TRUE) + 
  stat_poly_eq(aes(label = paste(
    # after_stat(eq.label),
    # after_stat(rr.label), 
    after_stat(p.value.label),
    sep = "~~~")),
    formula = y ~ x) + 
  geom_smooth(method = "lm", se = TRUE, color = "#886A9E", fill = "#f4e6ff", alpha = 0.5) + 
  geom_point() + 
  coord_cartesian(xlim = c(2018.9, 2022.1), ylim = c(0, 35), expand = 0) + 
  labs(title = "Lineage 1.3", x = "Sampling year", y = "Sampling numbers") + 
  theme_classic(base_size = 14)

ggsave("./plot/overview_plot/Lineage_1.3_sample_nums_by_year.pdf", width = 3, height = 3)

ggplot(data = dplyr::filter(summarize_data, lineage_level_2 == "1.2"), aes(x = collection_year, y = sample_nums, weight = yearly_total_samples)) + 
  geom_smooth(method = "lm", se = TRUE) + 
  stat_poly_eq(aes(label = paste(
    after_stat(p.value.label),
    sep = "~~~")),
    formula = y ~ x) + 
  geom_smooth(method = "lm", se = TRUE, color = "#9DD0CF", fill = "#e6ffff", alpha = 0.5) + 
  geom_point() + 
  coord_cartesian(xlim = c(2018.9, 2022.1), ylim = c(0, 35), expand = 0) + 
  labs(title = "Lineage 1.2", x = "Sampling year", y = "Sampling numbers") + 
  theme_classic(base_size = 14)

ggsave("./plot/overview_plot/Lineage_1.2_sample_nums_by_year.pdf", width = 3, height = 3)

ggplot(data = dplyr::filter(summarize_data, lineage_level_2 == "1.1"), aes(x = collection_year, y = sample_nums, weight = yearly_total_samples)) + 
  geom_smooth(method = "lm", se = TRUE) + 
  stat_poly_eq(aes(label = paste(
    after_stat(p.value.label),
    sep = "~~~")),
    formula = y ~ x) + 
  geom_smooth(method = "lm", se = TRUE, color = "#7EBB8F", fill = "#e6ffed", alpha = 0.5) + 
  labs(title = "Lineage 1.1", x = "Sampling year", y = "Sampling numbers") + 
  geom_point() + 
  theme_classic(base_size = 14)

ggsave("./plot/overview_plot/Lineage_1.1_sample_nums_by_year.pdf", width = 3, height = 3)
