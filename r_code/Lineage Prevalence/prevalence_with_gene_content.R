library(tidyverse)
library(broom)

setwd('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



prevalence_rate_by_year <- metadata.sample %>% 
  group_by(collection_year, lineage_level_2) %>% 
  summarise(covered_counties = n_distinct(county_gb), .groups = "drop")  %>%
  mutate(prevalence_rate = covered_counties / length(unique(metadata.sample$county_gb)))



## 使用以下方法还可分析基因家族拷贝数变异

# pangenes.panaroo <- read.csv('./data/panaroo_result/gene_presence_absence_roary.csv', header = TRUE, stringsAsFactors = FALSE)
# 
# datmat.counts <- column_to_rownames(pangenes.panaroo, var = "Gene")
# datmat.counts <- datmat.counts[, 14:ncol(datmat.counts)]
# datmat.counts <- as.data.frame(t(datmat.counts))
# rownames(datmat.counts) <- gsub("^X", "", rownames(datmat.counts))
# # datmat.counts <- datmat.counts[grepl("^20", rownames(datmat.counts)), ]
# datmat.counts <- as.data.frame(
#   apply(datmat.counts, c(1, 2), function(x) {
#     if (is.na(x) || x == "") {
#       return(0)  # NA或空字符串视为0
#     } else {
#       return(length(strsplit(x, ";")[[1]]))
#     }
#   })
# )
# datmat.counts <- datmat.counts[, colMeans(datmat.counts) != 0]
# 
# # same as datmat.panaroo totaly
# datmat.onehot <- datmat.counts
# datmat.onehot[datmat.onehot > 0] <- 1


normalize_datamat <- function(datmat) {
  colnames(datmat) <- datmat[1, ]
  datmat <- datmat[-1, ]

  rownames_datamat <- datmat$Gene
  datmat$Gene <- NULL

  datmat <- datmat %>% reframe(across(everything(), as.integer))
  rownames(datmat) <- rownames_datamat
  datmat <- as.data.frame(t(datmat))

  return(datmat)
}

datmat.onehot <- read.table('./data/panaroo_result_sensitive/gene_presence_absence.Rtab', header = FALSE, stringsAsFactors = FALSE)
datmat.onehot <- normalize_datamat(datmat.onehot)

softgene_names <- colnames(datmat.onehot)[colMeans(datmat.onehot) > 0.15 & colMeans(datmat.onehot) < 0.95]
raregene_names <- colnames(datmat.onehot)[colMeans(datmat.onehot) <= 0.15]



gene_content_data <- data.frame(
  assembly_accession = rownames(datmat.onehot), 
  family_content = rowSums(datmat.onehot), 
  soft_family_content  = rowSums(datmat.onehot[, softgene_names]), 
  rare_family_content  = rowSums(datmat.onehot[, raregene_names])
) %>% 
  right_join(., metadata.sample[, c('collection_year', 'assembly_accession', 'lineage_level_2')], by = 'assembly_accession') %>% 
  mutate(
    soft_family_content = soft_family_content / family_content,
    rare_family_content = rare_family_content / family_content,
  )

gene_content_summarize <- gene_content_data %>%
  group_by(collection_year, lineage_level_2) %>% 
  summarise(
    family_content_mean = mean(family_content, na.rm = TRUE),
    family_content_se   = sd(family_content, na.rm = TRUE) / sqrt(n()),
    
    soft_family_content_mean = mean(soft_family_content, na.rm = TRUE),
    soft_family_content_se   = sd(soft_family_content, na.rm = TRUE) / sqrt(n()),
    
    rare_family_content_mean = mean(rare_family_content, na.rm = TRUE),
    rare_family_content_se   = sd(rare_family_content, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  na.omit()

gene_content_data <- gene_content_data %>% 
  left_join(., prevalence_rate_by_year, by = c('lineage_level_2', 'collection_year'))

gene_content_summarize <- gene_content_summarize %>% 
  left_join(., prevalence_rate_by_year, by = c('lineage_level_2', 'collection_year'))



################################################################################
## 稀有基因

# 拟合线性模型
fit <- lm(rare_family_content_mean ~ covered_counties, data = gene_content_summarize)
fit_summary <- broom::glance(fit)
fit_coef <- broom::tidy(fit)

# 提取统计指标
r_squared <- round(fit_summary$r.squared, 2)
beta <- fit_coef$estimate[2]
p_value <- fit_coef$p.value[2]
p_label <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 2)))
# label_text <- paste0("R² = ", r_squared, ", β = ", beta, ", ", p_label)
label_text <- bquote(R^2 == .(r_squared) ~ "," ~ .(p_label))

print(beta)
# [1] 0.0001871233

ggplot(gene_content_summarize, aes(x = covered_counties, y = rare_family_content_mean, color = lineage_level_2, fill = lineage_level_2)) +
  # 添加外部拟合的回归线
  geom_smooth(data = gene_content_data, 
              aes(x = covered_counties, y = rare_family_content), 
              method = "lm", color = NA, fill = "grey80", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = rare_family_content_mean - rare_family_content_se,
                    ymax = rare_family_content_mean + rare_family_content_se),
                width = 0.5, size = 0.75) + 
  geom_point(size = 3, shape = 21, stroke = 1) +
  # 添加外部拟合的回归线
  geom_smooth(data = gene_content_data, 
              aes(x = covered_counties, y = rare_family_content), 
              method = "lm", color = "grey50", fill = NA, inherit.aes = FALSE, linewidth = 1.2) + 
  # 添加显著性文字
  annotate("text", x = Inf, y = Inf, label = label_text,
           hjust = 1.1, vjust = 1.5, size = 4) +
  labs(
    x = "Covered counties",
    y = "Rare-gene content"
  ) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_color_manual(name = "Lineage", values = c("1.1" = "#5e8c6b", "1.2" = "#7da6a5", "1.3" = "#6e5680")) + 
  scale_fill_manual(name = "Lineage", values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E")) + 
  theme_classic(base_size = 14)

# ggsave('./plot/Lineage Prevalence/rare_gene_content.pdf', width = 4.2, height = 3)

write.csv(gene_content_summarize, "./r_plot_data/figure_3c~gene_content_summarize.csv", row.names = FALSE)
write.csv(gene_content_data, "./r_plot_data/figure_3c~gene_content_data.csv", row.names = FALSE)

################################################################################
## 软基因

fit <- lm(soft_family_content_mean ~ covered_counties, data = gene_content_summarize)
fit_summary <- broom::glance(fit)
fit_coef <- broom::tidy(fit)

r_squared <- round(fit_summary$r.squared, 2)
beta <- fit_coef$estimate[2]
p_value <- fit_coef$p.value[2]
p_label <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 2)))
label_text <- bquote(R^2 == .(r_squared) ~ "," ~ .(p_label))

print(beta)
# [1] -0.0001964368

ggplot(gene_content_summarize, aes(x = covered_counties, y = soft_family_content_mean, color = lineage_level_2, fill = lineage_level_2)) +
  # 添加外部拟合的回归线
  geom_smooth(data = gene_content_data, 
              aes(x = covered_counties, y = soft_family_content), 
              method = "lm", color = NA, fill = "grey80", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = soft_family_content_mean - rare_family_content_se,
                    ymax = soft_family_content_mean + rare_family_content_se),
                width = 0.5, size = 0.75) +
  geom_point(size = 3, shape = 21, stroke = 1) +
  # 添加外部拟合的回归线
  geom_smooth(data = gene_content_data, 
              aes(x = covered_counties, y = soft_family_content), 
              method = "lm", color = "grey50", fill = NA, inherit.aes = FALSE, linewidth = 1.2) +
  # 添加显著性文字
  annotate("text", x = Inf, y = Inf, label = label_text,
           hjust = 1.1, vjust = 1.5, size = 4) +
  labs(
    x = "Covered counties",
    y = "Soft-gene content"
  ) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_color_manual(name = "Lineage", values = c("1.1" = "#5e8c6b", "1.2" = "#7da6a5", "1.3" = "#6e5680")) + 
  scale_fill_manual(name = "Lineage", values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E")) + 
  theme_classic(base_size = 14)

# ggsave('./plot/Lineage Prevalence/soft_gene_content.pdf', width = 4.2, height = 3)

write.csv(gene_content_summarize, "./r_plot_data/figure_3b~gene_content_summarize.csv", row.names = FALSE)
write.csv(gene_content_data, "./r_plot_data/figure_3b~gene_content_data.csv", row.names = FALSE)

################################################################################
## 总基因

fit <- lm(family_content_mean ~ covered_counties, data = gene_content_summarize)
fit_summary <- broom::glance(fit)
fit_coef <- broom::tidy(fit)

r_squared <- round(fit_summary$r.squared, 2)
beta <- fit_coef$estimate[2]
p_value <- fit_coef$p.value[2]
p_label <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 2)))
label_text <- bquote(R^2 == .(r_squared) ~ "," ~ .(p_label))

print(beta)
# [1] 0.02421601

ggplot(gene_content_summarize, aes(x = covered_counties, y = family_content_mean, color = lineage_level_2, fill = lineage_level_2)) +
  # 添加外部拟合的回归线
  geom_smooth(data = gene_content_data, 
              aes(x = covered_counties, y = family_content), 
              method = "lm", color = NA, fill = "grey80", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = family_content_mean - family_content_se,
                    ymax = family_content_mean + family_content_se),
                width = 0.5, size = 0.75) +
  geom_point(size = 3, shape = 21, stroke = 1) +
  # 添加外部拟合的回归线
  geom_smooth(data = gene_content_data, 
              aes(x = covered_counties, y = family_content), 
              method = "lm", color = "grey50", fill = NA, inherit.aes = FALSE, linewidth = 1.2) +
  # 添加显著性文字
  annotate("text", x = Inf, y = Inf, label = label_text,
           hjust = 1.1, vjust = 1.5, size = 4) +
  labs(
    x = "Covered Counties",
    y = "Genome Size"
  ) + 
  # scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_color_manual(name = "Lineage", values = c("1.1" = "#5e8c6b", "1.2" = "#7da6a5", "1.3" = "#6e5680")) + 
  scale_fill_manual(name = "Lineage", values = c("1.1" = "#7EBB8F", "1.2" = "#9DD0CF", "1.3" = "#886A9E")) + 
  theme_classic(base_size = 14)

# ggsave('./plot/Lineage Prevalence/genome_size.pdf', width = 4.2, height = 3)

write.csv(gene_content_summarize, "./r_plot_data/figure_3a~gene_content_summarize.csv", row.names = FALSE)
write.csv(gene_content_data, "./r_plot_data/figure_3a~gene_content_data.csv", row.names = FALSE)
