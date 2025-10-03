library(sf)        # 处理矢量数据
library(spdep)     # poly2nb
library(INLA)
library(tidyverse)
library(mgcv)      # gamm
library(nlme)      # 混合效应

setwd('/Users/xiyangli/Lab/Project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly')

load('./r_image/metadata_v3.RData')



sf_china <- read_sf('./source/中国_省.geojson') %>% dplyr::filter(name != '境界线')
## 修复无效几何
sf_china <- sf::st_make_valid(sf_china)
## 检查几何有效性
# sf::st_is_valid(sf_china)

## 定义空间邻接结构
# province_nb <- poly2nb(sf_china)
# nb2INLA("./source/province.adj", province_nb)
## 读取空间邻接结构
spatial_graph <- inla.read.graph(filename = "./source/province.adj")

## 设置省份 ID (和 graph 对应)
sf_china <- sf_china %>% mutate(ID_province = seq(1:nrow(.)))



national_data <- read.csv("./source/phsciencedata/布病.csv", header = TRUE)
# str(national_data)

colnames(national_data) <- c("short_name", "cases", "deaths", "incidence", "mortality", "year", "month")
national_data$short_name <- gsub(" ", "", national_data$short_name)
national_data$short_name <- gsub("省", "", national_data$short_name)
national_data$short_name <- gsub("市", "", national_data$short_name)
national_data <- national_data[national_data$short_name != "全国", ]
national_data <- national_data %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

# table(national_data$short_name)



mapping_data <- data.frame(
  full_name = c(
    "河南省",         "浙江省", "海南省", "台湾省",         "甘肃省",     "湖北省",         "天津市",   "内蒙古自治区",
    "广西壮族自治区", "江苏省", "山东省", "北京市",         "西藏自治区", "青海省",         "重庆市",   "澳门特别行政区",
    "吉林省",         "福建省", "广东省", "上海市",         "江西省",     "香港特别行政区", "黑龙江省", "山西省",
    "云南省",         "四川省", "安徽省", "宁夏回族自治区", "辽宁省",     "湖南省",         "陕西省",   "河北省",
    "新疆维吾尔自治区","贵州省"
  ),
  short_name = c(
    "河南", "浙江", "海南", "台湾", "甘肃", "湖北", "天津",   "内蒙古",
    "广西", "江苏", "山东", "北京", "西藏", "青海", "重庆",   "澳门",
    "吉林", "福建", "广东", "上海", "江西", "香港", "黑龙江", "山西",
    "云南", "四川", "安徽", "宁夏", "辽宁", "湖南", "陕西",   "河北",
    "新疆", "贵州"
  ),
  abbr_name = c(
    "HEN", "ZHJ", "HAN", NA,    "GAS", "HUB", "TAJ", "NEG",
    "GUX", "JSU", "SHD", "BEJ", "TIB", "QIH", "CHQ", NA,
    "JIL", "FUJ", "GUD", "SHH", "JXI", NA,    "HLL", "SHX",
    "YUN", "SCH", "ANH", "NXA", "LIA", "HUN", "SHA", "HEB",
    "XIN", "GUI"
  ),
  stringsAsFactors = FALSE
) %>% 
  dplyr::filter(!short_name %in% c("香港", "澳门", "台湾"))

national_data <- left_join(
  national_data, 
  mapping_data, 
  by = c("short_name" = "short_name")
)



data_fit <- left_join(
  national_data, 
  sf_china %>% st_drop_geometry() %>% select(name, ID_province), 
  by = c("full_name" = "name")
)

data_fit$year <- as.factor(data_fit$year)
data_fit$full_name <- as.factor(data_fit$full_name)



all_years <- c(2002, 2005, 2007, 2010, 2012, 2015, 2017, 2020)
result_list <- list()

for (yr in all_years) {
  # 读取数据
  file_path <- paste0("./source/food trade flows/A dataset of interprovincial food trade flows in China/",
                      yr, "matrix.xlsx")
  flow_data <- readxl::read_xlsx(file_path, sheet = "Beef and Mutton")
  flow_data <- column_to_rownames(flow_data, var = "Row")
  
  # 转换为矩阵（避免 rowSums/colSums 报错）
  flow_mat <- as.matrix(flow_data)
  
  bm_inner <- c()
  for (i in 1:nrow(flow_mat)) {
    bm_inner <- c(bm_inner, flow_data[i, i])
  }
  
  # 计算 rowSums 和 colSums
  bm_product <- rowSums(flow_mat, na.rm = TRUE)
  bm_output  <- rowSums(flow_mat, na.rm = TRUE) - bm_inner
  bm_input   <- colSums(flow_mat, na.rm = TRUE) - bm_inner
  
  # 合并为数据框（假设行名和列名一致，都是省份）
  df_out <- tibble(
    abbr_name  = names(bm_output),
    bm_product = as.numeric(bm_product), 
    bm_inner   = as.numeric(bm_inner), 
    bm_output  = as.numeric(bm_output),
    bm_input   = as.numeric(bm_input),
    year       = yr
  )
  
  result_list[[as.character(yr)]] <- df_out
}

# 合并所有年份
food_data <- bind_rows(result_list)

# 先确保 year 是数值型
food_data <- food_data %>% 
  mutate(year = as.integer(year))

# 补全年份并线性插值
food_data <- food_data %>%
  group_by(abbr_name) %>%
  tidyr::complete(year = 2002:2020) %>%   # 补全 2002-2020
  arrange(abbr_name, year) %>%
  mutate(
    bm_product = zoo::na.approx(bm_product, year, na.rm = FALSE),
    bm_inner   = zoo::na.approx(bm_inner, year, na.rm = FALSE),
    bm_output  = zoo::na.approx(bm_output, year, na.rm = FALSE),
    bm_input   = zoo::na.approx(bm_input,  year, na.rm = FALSE)
  ) %>%
  ungroup()

food_data$year <- as.factor(food_data$year)

data_fit <- left_join(
  data_fit, 
  food_data, 
  by = c("abbr_name" = "abbr_name", "year" = "year")
)










mod_gamm_0 <- mgcv::bam(
  formula = incidence ~ s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"), 
  data = data_fit,
  family = quasipoisson()
)

summary(mod_gamm_0)$r.sq

mod_gamm_1 <- mgcv::bam(
  formula = incidence ~ s(bm_product, k = 5) + s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"),
  data = data_fit,
  family = quasipoisson()
)

summary(mod_gamm_1)$r.sq
anova(mod_gamm_0, mod_gamm_1)

# plot(mod_gamm_1)


mod_gamm_2 <- mgcv::bam(
  formula = incidence ~ s(bm_product, k = 5) + s(bm_output, k = 5) + s(bm_input, k = 5) + s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"),
  data = data_fit,
  family = quasipoisson()
)

summary(mod_gamm_2)$r.sq
# plot(mod_gamm_2)

mod_gamm_3 <- mgcv::bam(
  formula = incidence ~ s(bm_product, k = 5) + s(bm_output, k = 5) + s(bm_input, k = 5, by = full_name) + s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"),
  data = data_fit,
  family = quasipoisson()
)

summary(mod_gamm_3)$r.sq
anova(mod_gamm_3, mod_gamm_2)
# plot(mod_gamm_3)



data_pred <- data_fit
data_pred$pred_0 <- predict(mod_gamm_0, newdata = data_fit, type = "response")
data_pred$pred_1 <- predict(mod_gamm_1, newdata = data_fit, type = "response")
data_pred$pred_2 <- predict(mod_gamm_2, newdata = data_fit, type = "response")
data_pred$pred_3 <- predict(mod_gamm_3, newdata = data_fit, type = "response")

data_pred %>%
  # dplyr::filter(short_name == "江苏") %>%
  dplyr::filter(short_name == "云南") %>%
  ggplot(., aes(x = date)) +
  geom_line(aes(y = incidence), color = "black") + 
  geom_line(aes(y = pred_0), color = "#c7ecee", linewidth = 1) + 
  geom_line(aes(y = pred_1), color = "#7ed6df", linewidth = 1) + 
  geom_line(aes(y = pred_2), color = "#1289a7", linewidth = 1) + 
  geom_line(aes(y = pred_3), color = "#3c6382", linewidth = 1) + 
  theme_minimal()



library(marginaleffects)

eff_yunnan_input <- predictions(
  model = mod_gamm_3,
  variables = list(bm_input = seq(min(data_fit[data_fit$full_name == "云南省", ]$bm_input),
                                  max(data_fit[data_fit$full_name == "云南省", ]$bm_input),
                                  length.out = 100)),
  newdata = datagrid(full_name = "云南省") # 固定云南省，其余自动边际化
)

eff_yunnan_product <- predictions(
  model = mod_gamm_3,
  variables = list(bm_product = seq(min(data_fit[data_fit$full_name == "云南省", ]$bm_product), 
                                    max(data_fit[data_fit$full_name == "云南省", ]$bm_product), 
                                    length.out = 100)),
  newdata = datagrid(full_name = "云南省") 
)

eff_jiangsu_input <- predictions(
  model = mod_gamm_3,
  variables = list(bm_input = seq(min(data_fit[data_fit$full_name == "江苏省", ]$bm_input),
                                  max(data_fit[data_fit$full_name == "江苏省", ]$bm_input),
                                  length.out = 100)),
  newdata = datagrid(full_name = "江苏省") 
)

eff_jiangsu_product <- predictions(
  model = mod_gamm_3,
  variables = list(bm_product = seq(min(data_fit[data_fit$full_name == "江苏省", ]$bm_product), 
                                    max(data_fit[data_fit$full_name == "江苏省", ]$bm_product), 
                                    length.out = 100)),
  newdata = datagrid(full_name = "江苏省") 
)

ggplot(eff_jiangsu_input, aes(x = bm_input, y = estimate)) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()

ggplot(eff_yunnan_input, aes(x = bm_input, y = estimate)) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()


ggplot(eff_jiangsu_product, aes(x = bm_product, y = estimate)) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()

ggplot(eff_yunnan_product, aes(x = bm_product, y = estimate)) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()

















mod_gamm_g0 <- mgcv::bam(
  formula = log(incidence) ~ s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"), 
  data = data_fit, 
  family = gaussian()
)

summary(mod_gamm_g0)$r.sq

mod_gamm_g1 <- mgcv::bam(
  formula = log(incidence) ~ s(bm_product) + s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"), 
  data = data_fit, 
  family = gaussian()
)

summary(mod_gamm_g1)$r.sq

mod_gamm_g2 <- mgcv::bam(
  formula = log(incidence) ~ s(bm_product) + s(bm_output) + s(bm_input) + s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"), 
  data = data_fit, 
  family = gaussian()
)

summary(mod_gamm_g2)$r.sq



mod_gamm_g3 <- mgcv::bam(
  formula = log(incidence) ~ s(bm_product, k = 5) + s(bm_input, by = full_name) + s(bm_output, by = full_name) + s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"), 
  data = data_fit, 
  family = gaussian()
)

summary(mod_gamm_g3)$r.sq

# mod_gamm_g3.1 <- mgcv::bam(
#   formula = log(incidence) ~ s(bm_product) + s(bm_input, by = full_name) + s(full_name, bs = "re") + s(year, bs = "re") + s(month, bs = "cr"), 
#   data = data_fit, 
#   family = gaussian()
# )
# 
# summary(mod_gamm_g3.1)$r.sq
# 
# anova(mod_gamm_g3.1, mod_gamm_g3)

# par(mfrow = c(2, 3))
# plot(mod_gamm_2, shade = TRUE)
# 
# par(mfrow = c(2, 2))
# mgcv::gam.check(mod_gamm$gam) 



data_pred <- data_fit
data_pred$pred_g0 <- exp(predict(mod_gamm_g0, newdata = data_fit, type = "response"))
data_pred$pred_g1 <- exp(predict(mod_gamm_g1, newdata = data_fit, type = "response"))
data_pred$pred_g2 <- exp(predict(mod_gamm_g2, newdata = data_fit, type = "response"))
data_pred$pred_g3 <- exp(predict(mod_gamm_g3, newdata = data_fit, type = "response"))

data_pred %>% 
  # dplyr::filter(short_name == "江苏") %>%
  dplyr::filter(short_name == "云南") %>%
  ggplot(., aes(x = date)) +
  geom_line(aes(y = incidence), color = "black", linewidth = 1) + 
  # geom_line(aes(y = pred_g0), color = "#c7ecee", linewidth = 1) +
  geom_line(aes(y = pred_g1), color = "#ce8892", linewidth = 0.75) +
  # geom_line(aes(y = pred_g2), color = "#1289a7", linewidth = 1) +
  geom_line(aes(y = pred_g3), color = "#7ed6df", linewidth = 0.75) + 
  theme_classic()

data_pred %>% 
  # dplyr::filter(short_name == "江苏") %>%
  dplyr::filter(short_name == "云南") %>%
  ggplot(., aes(x = date)) +
  geom_line(aes(y = bm_input), color = "#3c6382", linewidth = 1) + 
  theme_minimal()

library(marginaleffects)

eff_yunnan_input <- predictions(
  model = mod_gamm_g3,
  variables = list(bm_input = seq(min(data_fit[data_fit$full_name == "云南省", ]$bm_input),
                                  max(data_fit[data_fit$full_name == "云南省", ]$bm_input),
                                  length.out = 100)),
  newdata = datagrid(full_name = "云南省") # 固定云南省，其余自动边际化
)

eff_yunnan_product <- predictions(
  model = mod_gamm_g3,
  variables = list(bm_product = seq(min(data_fit[data_fit$full_name == "云南省", ]$bm_product), 
                                    max(data_fit[data_fit$full_name == "云南省", ]$bm_product), 
                                    length.out = 100)),
  newdata = datagrid(full_name = "云南省") 
)

eff_yunnan_output <- predictions(
  model = mod_gamm_g3,
  variables = list(bm_output = seq(min(data_fit[data_fit$full_name == "云南省", ]$bm_output), 
                                   max(data_fit[data_fit$full_name == "云南省", ]$bm_output), 
                                   length.out = 100)),
  newdata = datagrid(full_name = "云南省") 
)

eff_jiangsu_input <- predictions(
  model = mod_gamm_g3,
  variables = list(bm_input = seq(min(data_fit[data_fit$full_name == "江苏省", ]$bm_input),
                                  max(data_fit[data_fit$full_name == "江苏省", ]$bm_input),
                                  length.out = 100)),
  newdata = datagrid(full_name = "江苏省") 
)

eff_jiangsu_product <- predictions(
  model = mod_gamm_g3,
  variables = list(bm_product = seq(min(data_fit[data_fit$full_name == "江苏省", ]$bm_product), 
                                    max(data_fit[data_fit$full_name == "江苏省", ]$bm_product), 
                                    length.out = 100)),
  newdata = datagrid(full_name = "江苏省") 
)

eff_jiangsu_output <- predictions(
  model = mod_gamm_g3,
  variables = list(bm_output = seq(min(data_fit[data_fit$full_name == "江苏省", ]$bm_output), 
                                   max(data_fit[data_fit$full_name == "江苏省", ]$bm_output), 
                                   length.out = 100)),
  newdata = datagrid(full_name = "江苏省") 
)

ggplot(eff_jiangsu_input, aes(x = bm_input, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()

ggplot(eff_yunnan_input, aes(x = bm_input, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()

ggplot(eff_jiangsu_output, aes(x = bm_output, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()

ggplot(eff_yunnan_output, aes(x = bm_output, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) +
  theme_classic()

ggplot(eff_jiangsu_product, aes(x = bm_product, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()

ggplot(eff_yunnan_product, aes(x = bm_product, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
  theme_classic()









# 从原始数据里提取云南省 2004-2020 的年份与 bm_input 对应值
dat_yunnan <- data_fit %>%
  filter(full_name == "云南省") %>%
  select(year, bm_input) %>%
  unique.data.frame()

# 构造 newdata：固定云南省，并带入逐年的 bm_input
newdat <- datagrid(
  model = mod_gamm_g3, 
  year = dat_yunnan$year,
  # month = seq(1, 12, 1), 
  full_name = "云南省"
)

newdat$bm_input <- NULL
newdat <- left_join(newdat, dat_yunnan, by = "year")

# 计算预测值（边际效应）
eff_yunnan_input_year <- predictions(
  model = mod_gamm_g3,
  newdata = newdat
)

# eff_yunnan_input_year %>% 
#   group_by(year) %>% 
#   summarise(estimate = mean(estimate), 
#             conf.low = mean(conf.low), 
#             conf.high = mean(conf.high), 
#             .groups = "drop") %>% 
#   ggplot(., aes(x = as.numeric(year), y = exp(estimate))) +
#   geom_line(color = "#0a3d62") +
#   geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
#   theme_classic()

ggplot(eff_yunnan_input_year, aes(x = as.numeric(year) + 2003, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) +
  theme_classic()



dat_yunnan_product <- data_fit %>%
  filter(full_name == "云南省") %>%
  select(year, bm_product) %>%
  unique.data.frame()

# 构造 newdata：固定云南省，并带入逐年的 bm_input
newdat_product <- datagrid(
  model = mod_gamm_g3, 
  year = dat_yunnan$year,
  full_name = "云南省"
)

newdat_product$bm_product <- NULL
newdat_product <- left_join(newdat_product, dat_yunnan_product, by = "year")

# 计算预测值（边际效应）
eff_yunnan_product_year <- predictions(
  model = mod_gamm_g3,
  newdata = newdat_product
)

ggplot(eff_yunnan_product_year, aes(x = as.numeric(year) + 2003, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) +
  theme_classic()



ggplot() +
  geom_line(data = eff_yunnan_product_year, 
            aes(x = as.numeric(year) + 2003, y = exp(estimate)), 
            color = "#0a3d62") +
  geom_ribbon(data = eff_yunnan_product_year,
              aes(x = as.numeric(year) + 2003, ymin = exp(conf.low), ymax = exp(conf.high)),
              fill = "#778ca3", alpha = 0.2) +
  geom_line(data = eff_yunnan_input_year, 
            aes(x = as.numeric(year) + 2003, y = exp(estimate)), 
            color = "#0a3d62") +
  geom_ribbon(data = eff_yunnan_input_year, 
              aes(x = as.numeric(year) + 2003, ymin = exp(conf.low), ymax = exp(conf.high)), 
              fill = "#778ca3", alpha = 0.2) +
  theme_classic()





# save.image('./r_image/spat_model.RData')
load('./r_image/spat_model.RData')







p1 <- data_pred %>% 
  dplyr::filter(short_name == "云南") %>%
  ggplot(., aes(x = date)) +
  geom_line(aes(y = incidence), color = "grey80", linewidth = 1) + 
  geom_line(aes(y = pred_g1), color = "#546de5", linewidth = 1) + 
  theme_classic(base_size = 14) + 
  labs(x = "Date", y = "Incidence (1/100,000)")

p2 <- data_pred %>% 
  dplyr::filter(short_name == "云南") %>%
  ggplot(., aes(x = date)) +
  geom_line(aes(y = incidence), color = "grey80", linewidth = 1) + 
  geom_line(aes(y = pred_g3), color = "#e55039", linewidth = 1) + 
  theme_classic(base_size = 14) + 
  labs(x = "Date", y = "Incidence (1/100,000)")

# data_pred %>% 
#   dplyr::filter(short_name == "云南") %>%
#   ggplot(., aes(x = date)) +
#   geom_abline(intercept = 0, linetype = 2) + 
#   geom_point(aes(x = incidence, y = pred_g1), color = "#546de5", size = 3, alpha = 0.5) + 
#   theme_classic(base_size = 14)
# 
# data_pred %>% 
#   dplyr::filter(short_name == "云南") %>%
#   ggplot(., aes(x = date)) +
#   geom_abline(intercept = 0, linetype = 2) + 
#   geom_point(aes(x = incidence, y = pred_g3), color = "#e55039", size = 3, alpha = 0.5) + 
#   theme_classic(base_size = 14)

p3 <- ggplot(eff_yunnan_input, aes(x = bm_input, y = exp(estimate))) +
  geom_line(color = "#0a3d62") +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
  theme_classic(base_size = 14) + 
  labs(x = "Total inflow", y = "Partial effect")

p4 <- ggplot(eff_yunnan_input_year, aes(x = as.numeric(year) + 2003, y = exp(estimate))) +
  geom_line(color = "#0a3d62") + 
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#778ca3", alpha = 0.2) + 
  geom_vline(xintercept = 2015, linetype = 2) + 
  theme_classic(base_size = 14) + 
  labs(x = "Date", y = "Marginal effect of total inflow")



load('./r_image/metadata_v3.RData')

profession_prop <- metadata.cases %>% 
  mutate(profession = case_when(
    profession == "农民"         ~ "1",
    profession == "牧民"         ~ "1",
    profession == "散居儿童"     ~ "3",
    profession == "教师"         ~ "3",
    profession == "学生"         ~ "3",
    profession == "医务人员"     ~ "3",
    profession == "商业服务"     ~ "3",
    profession == "工人"         ~ "3",
    profession == "幼托儿童"     ~ "3",
    profession == "民工"         ~ "3",
    profession == "家务及待业"   ~ "3",
    profession == "离退人员"     ~ "3",
    profession == "干部职员"     ~ "3",
    profession == "餐饮食品业"   ~ "1",
    profession == "其他"         ~ "3",
    profession == "不详"         ~ "3",
    TRUE                         ~ NA_character_
  )) %>%
  mutate(
    year = year(onset_date)
  ) %>%
  # dplyr::filter(year >= 2010) %>% 
  group_by(profession, year) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  na.omit() %>%
  group_by(year) %>% 
  mutate(total = sum(count), 
         prop  = count / total) %>% 
  dplyr::filter(profession == "3") %>%
  ungroup()

profession_prop <- add_row(profession_prop, data.frame(profession = "3", year = 2011, count = 0, total = 16, prop = 0))






ggplot() + 
  # geom_smooth(method = 'lm', formula = 'y ~ x', color = "grey10", fill = "#e7e7eb") +
  geom_smooth(data = profession_prop,
              aes(x = year, y = prop, weight = total),
              method = 'gam', formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5), color = "#e55039", fill = "#f8c291", alpha = 0.2) +
  # geom_smooth(data = profession_prop, 
  #             aes(x = year, y = prop, weight = total), 
  #             method = 'loess', formula = y ~ x, color = "#0a3d62", fill = "#778ca3") + 
  # geom_point(data = profession_prop, 
  #            aes(x = year, y = prop, size = count), 
  #            color = "#0a3d62", alpha = 0.8) + 
  geom_vline(xintercept = 2015, linetype = 2) + 
  coord_cartesian(ylim = c(0, 0.3), x = c(2008, 2020)) + 
  theme_classic(base_size = 14) + 
  geom_line(data = eff_yunnan_input_year, 
            aes(x = as.numeric(year) + 2003, y = exp(estimate)), 
            color = "#0a3d62", linewidth = 1) +
  geom_ribbon(data = eff_yunnan_input_year, 
              aes(x = as.numeric(year) + 2003, ymin = exp(conf.low), ymax = exp(conf.high)), 
              fill = "#778ca3", alpha = 0.2)


BETA <- 0.8
ALPHA <- 0.12
p7 <- ggplot() + 
  # 副轴：畜牧贸易边际效应（缩放后）
  geom_ribbon(
    data = eff_yunnan_input_year, 
    aes(x = as.numeric(year) + 2003, 
        ymin = exp(conf.low) * BETA + ALPHA, 
        ymax = exp(conf.high) * BETA + ALPHA), 
    fill = "#778ca3", alpha = 0.2
  ) +
  
  # 主轴：职业比例
  geom_smooth(
    data = profession_prop,
    aes(x = year, y = prop, weight = total),
    method = 'gam',
    formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5),
    color = "#e55039",
    fill = "#f8c291",
    alpha = 0.2
  ) + 
  geom_vline(xintercept = 2015, linetype = 2) + 
  
  # 副轴：畜牧贸易边际效应（缩放后）
  geom_line(
    data = eff_yunnan_input_year, 
    aes(x = as.numeric(year) + 2003, y = exp(estimate) * BETA + ALPHA), 
    color = "#0a3d62", linewidth = 1
  ) +
  
  # 坐标轴设置
  scale_y_continuous(
    name = "Annual proportion of cases in \nnon-agro-pastoral occupations",
    sec.axis = sec_axis(
      ~ (. - ALPHA) / BETA, 
      name = "Marginal effect of livestock \ntrade volume on incidence"
    )
  ) +
  scale_x_continuous(
    name = "Date",
    breaks = seq(2008, 2020, 2) # 只显示整数
  ) +
  coord_cartesian(ylim = c(0.12, 0.35), xlim = c(2008, 2020), expand = 0) +
  theme_classic(base_size = 14)






metadata.sample %>%
  mutate(profession = case_when(
    profession == "农民"         ~ "a) Farmer",
    profession == "牧民"         ~ "a) Herder",
    profession == "散居儿童"     ~ "c) Non-institutional Child",
    profession == "教师"         ~ "c) Teacher",
    profession == "学生"         ~ "c) Student",
    profession == "医务人员"     ~ "c) Doctor",
    profession == "商业服务"     ~ "c) Business Services",
    profession == "工人"         ~ "c) Worker",
    profession == "幼托儿童"     ~ "c) Children in childcare",
    profession == "民工"         ~ "c) Migrant workers",
    profession == "家务及待业"   ~ "c) Homemaker / Unemployed",
    profession == "离退人员"     ~ "c) Retired personnel",
    profession == "干部职员"     ~ "c) Office staff",
    profession == "餐饮食品业"   ~ "c) Catering / Food",
    profession == "其他"         ~ "b) Unknown / Other",
    profession == "不详"         ~ "b) Unknown / Other",
    TRUE                         ~ NA_character_
  )) %>%
  mutate(
    year = year(onset_date)
  ) %>%
  group_by(profession, year) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  group_by(year) %>% 
  mutate(total = sum(count), .groups = "drop", 
         prop  = count / total) %>% 
  ggplot(., aes(x = year, y = prop, fill = profession)) + 
  geom_col(position = "stack", color = "black") + 
  scale_fill_manual(
    values = c(
      "a) Farmer"                   = "#34ace0",
      "a) Herder"                   = "#63cdda",
      "c) Business Services"        = "#f8e58c",
      "c) Catering / Food"          = "#fcd575",
      "c) Children in childcare"    = "#fbca4d",
      "c) Doctor"                   = "#fcc800",
      "c) Homemaker / Unemployed"   = "#fabf14",
      "c) Migrant workers"          = "#f6ad49",
      "c) Non-institutional Child"  = "#f39800",
      "c) Office staff"             = "#f08300",
      "c) Retired personnel"        = "#ee7800",
      "c) Student"                  = "#ec6800",
      "c) Teacher"                  = "#ea5506",
      "c) Worker"                   = "#d3381c",
      "b) Unknown / Other"          = "grey80"
    )
  ) + 
  labs(x = "Year", y = "Occupational proportion of sampled cases") + 
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14)

p5 <- metadata.sample %>%
  mutate(profession = case_when(
    profession == "农民"         ~ "a) Farmer",
    profession == "牧民"         ~ "a) Herder",
    profession == "散居儿童"     ~ "c) Non-institutional Child",
    profession == "教师"         ~ "c) Teacher",
    profession == "学生"         ~ "c) Student",
    profession == "医务人员"     ~ "c) Doctor",
    profession == "商业服务"     ~ "c) Business Services",
    profession == "工人"         ~ "c) Worker",
    profession == "幼托儿童"     ~ "c) Children in childcare",
    profession == "民工"         ~ "c) Migrant workers",
    profession == "家务及待业"   ~ "c) Homemaker / Unemployed",
    profession == "离退人员"     ~ "c) Retired personnel",
    profession == "干部职员"     ~ "c) Office staff",
    profession == "餐饮食品业"   ~ "c) Catering / Food",
    profession == "其他"         ~ "b) Unknown / Other",
    profession == "不详"         ~ "b) Unknown / Other",
    TRUE                         ~ NA_character_
  )) %>%
  group_by(profession, lineage_level_2) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(., aes(x = lineage_level_2, y = count, fill = profession)) + 
  geom_col(position = "stack", color = "black") + 
  scale_fill_manual(
    values = c(
      "a) Farmer"                   = "#34ace0",
      "a) Herder"                   = "#63cdda",
      "c) Business Services"        = "#f8e58c",
      "c) Catering / Food"          = "#fcd575",
      "c) Children in childcare"    = "#fbca4d",
      "c) Doctor"                   = "#fcc800",
      "c) Homemaker / Unemployed"   = "#fabf14",
      "c) Migrant workers"          = "#f6ad49",
      "c) Non-institutional Child"  = "#f39800",
      "c) Office staff"             = "#f08300",
      "c) Retired personnel"        = "#ee7800",
      "c) Student"                  = "#ec6800",
      "c) Teacher"                  = "#ea5506",
      "c) Worker"                   = "#d3381c",
      "b) Unknown / Other"          = "grey80"
    )
  ) + 
  labs(x = "Year", y = "Number of occupations in sampled cases") + 
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14)

p6 <- metadata.cases %>%
  mutate(profession = case_when(
    profession == "农民"         ~ "a) Farmer",
    profession == "牧民"         ~ "a) Herder",
    profession == "散居儿童"     ~ "c) Non-institutional Child",
    profession == "教师"         ~ "c) Teacher",
    profession == "学生"         ~ "c) Student",
    profession == "医务人员"     ~ "c) Doctor",
    profession == "商业服务"     ~ "c) Business Services",
    profession == "工人"         ~ "c) Worker",
    profession == "幼托儿童"     ~ "c) Children in childcare",
    profession == "民工"         ~ "c) Migrant workers",
    profession == "家务及待业"   ~ "c) Homemaker / Unemployed",
    profession == "离退人员"     ~ "c) Retired personnel",
    profession == "干部职员"     ~ "c) Office staff",
    profession == "餐饮食品业"   ~ "c) Catering / Food",
    profession == "其他"         ~ "b) Unknown / Other",
    profession == "不详"         ~ "b) Unknown / Other",
    TRUE                         ~ NA_character_
  )) %>%
  mutate(
    year = year(onset_date)
  ) %>%
  group_by(profession, year) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  # na.omit() %>% 
  group_by(year) %>% 
  mutate(total = sum(count), .groups = "drop", 
         prop  = count / total) %>% 
  ggplot(., aes(x = year, y = prop, fill = profession)) + 
  geom_col(position = "stack", color = "black") + 
  scale_fill_manual(
    values = c(
      "a) Farmer"                   = "#34ace0",
      "a) Herder"                   = "#63cdda",
      "c) Business Services"        = "#f8e58c",
      "c) Catering / Food"          = "#fcd575",
      "c) Children in childcare"    = "#fbca4d",
      "c) Doctor"                   = "#fcc800",
      "c) Homemaker / Unemployed"   = "#fabf14",
      "c) Migrant workers"          = "#f6ad49",
      "c) Non-institutional Child"  = "#f39800",
      "c) Office staff"             = "#f08300",
      "c) Retired personnel"        = "#ee7800",
      "c) Student"                  = "#ec6800",
      "c) Teacher"                  = "#ea5506",
      "c) Worker"                   = "#d3381c",
      "b) Unknown / Other"          = "grey80"
    )
  ) + 
  labs(x = "Year", y = "Occupational proportion of report cases") + 
  coord_cartesian(xlim = c(2007.5, NA), expand = FALSE) + 
  theme_classic(base_size = 14)

library(patchwork)

(p1 / p2 / (p3 | p4)) / (p5 | p6) / p7



(
  p1 / 
  p2 / 
  (p3 | p4) + plot_layout(widths = c(0.3, 0.7))
) / 
  (p5 | p6) / 
  p7 + plot_layout(guides = 'collect')

ggsave("./plot/trade_analysis/merge_plot.pdf", width = 12, height = 16)


