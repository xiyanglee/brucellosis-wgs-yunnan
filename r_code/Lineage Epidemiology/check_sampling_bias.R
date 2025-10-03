library(tidyverse)
library(cowplot)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



summarize_per_month <- metadata.cases %>%
  mutate(month = floor_date(onset_date, "month")) %>%  # 提取每个样本的月份（日期取整到月）
  group_by(month) %>%
  summarise(
    total_samples = n(),                   # 总样本数
    sequenced_samples = sum(!is.na(tip))     # 测序样本数（tip 非 NA）
  ) %>% 
  # dplyr::filter(sequenced_samples > 0) %>% 
  arrange(month)

m <- lm(sequenced_samples ~ total_samples - 1, data = summarize_per_month %>% filter(year(month) >= 2019))
summary_model <- summary(m)

# 提取 R² 和 p-value
r_squared <- summary_model$r.squared
p_value <- coef(summary_model)["total_samples", "Pr(>|t|)"]

# 绘制散点图并添加过原点的回归线
summarize_per_month %>% 
  filter(year(month) >= 2019) %>% 
  ggplot(aes(x = total_samples, y = sequenced_samples)) + 
  geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, color = "#3a5b52", linetype = 2) + 
  # geom_point(aes(fill = month), color = 'black', alpha = 1, size = 2, shape = 21) + # '#779649'
  geom_point(fill = '#779649', color = 'black', alpha = 1, size = 3, shape = 21) + 
  labs(x = "Total cases per month", y = "Sequenced samples per month") + 
  theme_classic(base_size = 14) + 
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) + 
  # theme(text = element_text(family = "Arial")) + 
  # scale_fill_gradientn(colors = c('#e0ebaf', '#316745'),
  #                      labels = function(x) format(as.POSIXct(x, origin = '1970-01-01'), '%Y')) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  # 在图中标注 R² 和 p-value（右上角）
  annotate("text", 
           x = 125, y = 3.5, 
           label = sprintf("R² = %.3f\np = %.3g", r_squared, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5)

# ggsave('./plot/check_sampling_bias/point~total_samples&sequenced_samples~per_month.pdf' , width = 4, height = 4)



stacked_seq <- metadata.cases %>%
  filter(!is.na(tip)) %>%  # 只保留有测序的样本
  mutate(month = floor_date(onset_date, "month")) %>%
  group_by(month, lineage_level_2) %>%
  summarise(sequenced_samples = n(), .groups = "drop")

p_inset <- summarize_per_month %>% 
  filter(year(month) >= 2019) %>% 
  ggplot(aes(x = total_samples, y = sequenced_samples)) + 
  geom_point(color = '#6b798e', size = 2, shape = 16, alpha = 0.5) + 
  geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, color = "#003460", linewidth = 0.75) + 
  labs(x = "Total cases per month", y = "Sequenced samples per month") + 
  theme_classic(base_size = 12) + 
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  # 在图中标注 R² 和 p-value（右上角）
  theme(
    panel.background = element_blank(),  # 面板背景透明
    plot.background = element_blank()    # 图整体背景透明
  ) + 
  annotate("text", 
           x = 125, y = 3.5, 
           label = sprintf("R² = %.3f\np = %.3g", r_squared, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5)

# 计算 total_samples 映射到左侧 y 轴的数值（变换使得范围匹配 0 ~ 6）
p_main <- summarize_per_month %>% 
  mutate(trans_total_samples = (1 / 10) * (summarize_per_month$total_samples + 0)) %>% 
  ggplot(aes(x = month)) +
  # 散点图：total_samples（经过线性变换）
  geom_line(aes(y = trans_total_samples), color = "#0E2350", linewidth = 0.5) + 
  geom_point(aes(y = trans_total_samples), color = "white", size = 3, shape = 16) +
  geom_point(aes(y = trans_total_samples), color = "#0E2350", fill = "white", size = 1, shape = 21, stroke = 1) +
  # 左侧 y 轴：sequenced_samples 柱状图
  # geom_bar(aes(y = sequenced_samples), stat = "identity", fill = "#8EA5D6", alpha = 1, width = 31 * 86400) + 
  geom_col(data = stacked_seq,
           aes(x = month, y = sequenced_samples, fill = lineage_level_2),
           width = 31 * 86400, alpha = 1, position = "stack") + 
  scale_fill_manual(values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E")) + 
  scale_y_continuous(
    name = "Sequenced Samples per month",
    expand = expansion(mult = c(0, 0.05)),
    # 右侧 y 轴：逆变换还原 total_samples 的刻度
    sec.axis = sec_axis(~ . * (10 / 1) - 0, name = "Total Samples per month"), 
  ) +
  labs(x = "Year") + 
  theme_bw(base_size = 14) +
  coord_cartesian(ylim = c(0, 14), xlim = c(as.POSIXct("2010-01-01"), as.POSIXct("2022-10-15")), expand = FALSE) +
  theme(
    axis.title.y = element_text(color = "#003460"),
    axis.title.y.right = element_text(color = "#0E2350"), 
    panel.grid = element_blank(), 
    legend.position = "none"
  )

ggdraw() +
  draw_plot(p_main) +
  draw_plot(p_inset, x = 0.13, y = 0.4, width = 0.45, height = 0.52)

# ggsave('./plot/check_sampling_bias/cowplot~year&sequenced_samples&total_samples~per_month.pdf' , width = 4.5, height = 4)
rm(m, summary_model, r_squared, p_value)

write.csv(summarize_per_month, "./r_plot_data/figure_2a~summarize_per_month.csv", row.names = FALSE)


summarize_per_city <- metadata.cases %>% 
  filter(year(onset_date) >= 2019) %>% 
  mutate(city = ifelse(str_detect(address, "省"),
                       # When "省" exists, extract the text between "省" and the first occurrence of "市" or "州"
                       stringr::str_extract(address, "(?<=省)[^市州]+[市州]"),
                       # Otherwise, assume the address starts with the city name ending with 市 or 州
                       stringr::str_extract(address, "^[^市州]+[市州]"))) %>%
  group_by(city) %>%
  summarise(total_samples = n(),
            sequenced_samples = sum(!is.na(tip))) %>% 
  # dplyr::filter(sequenced_samples > 0) %>% 
  ungroup()

m <- lm(sequenced_samples ~ total_samples - 1, data = summarize_per_city %>% filter(sequenced_samples > 0))
summary_model <- summary(m)

# 提取 R² 和 p-value
r_squared <- summary_model$r.squared
p_value <- coef(summary_model)["total_samples", "Pr(>|t|)"]

# 绘制散点图并添加过原点的回归线
ggplot(summarize_per_city %>% filter(sequenced_samples > 0), aes(x = log10(total_samples + 1), y = log10(sequenced_samples + 1))) + 
  geom_point(color = 'grey', alpha = 0.5, size = 2) + 
  geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, color = "blue") + 
  labs(x = "Total Samples", y = "Sequenced Samples") + 
  theme_classic() + 
  # coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) + 
  # 在图中标注 R² 和 p-value（右上角）
  annotate("text", 
           x = 3, y = 1, 
           label = sprintf("R² = %.3f\np = %.3g", r_squared, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5)

ggplot(summarize_per_city, aes(x = total_samples, y = sequenced_samples)) + 
  geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, color = "#3a5b52", linetype = 2) + 
  geom_point(fill = '#779649', color = 'black', alpha = 1, size = 3, shape = 21) + 
  labs(x = "Total cases per month", y = "Sequenced samples per month") + 
  theme_classic(base_size = 14) + 
  coord_cartesian(xlim = c(1, NA), ylim = c(1, NA)) +
  scale_x_log10(expand = expansion(mult = c(0, 0.05)),
                # breaks = scales::trans_breaks("log10", function(x) 10^x),
                breaks = c(1, 10, 100, 1000), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(expand = expansion(mult = c(0, 0.05)),
                # breaks = scales::trans_breaks("log10", function(x) 10^x),
                breaks = c(1, 10, 10^1.5), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotate("text", 
           x = 20, y = 10, 
           label = sprintf("R² = %.3f\np = %.3g", r_squared, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5) + 
  annotation_logticks()

# ggsave('./plot/check_sampling_bias/point~total_samples&sequenced_samples~per_city.pdf' , width = 4, height = 4)
rm(m, summary_model, r_squared, p_value)

