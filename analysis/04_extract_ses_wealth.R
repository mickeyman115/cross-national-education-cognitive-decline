##############################################################################
# 13. 提取财富、收入与童年 SES (13_extract_ses.R)
# 修正版本：处理 CHARLS 的宽转长拆分问题，提取完整童年 SES，排除 SHARE 2017
##############################################################################

library(dplyr)
library(haven)
library(tidyr)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
out_dir <- PROJECT_DIR
df_main <- readRDS(file.path(out_dir, "clean_longitudinal_data.rds"))

# 获取分析基线
baseline_waves <- df_main %>% filter(!is.na(cogtot)) %>% 
  group_by(id) %>% 
  summarize(base_wave = min(wave)) %>%
  ungroup()

cat("================================================================\n")
cat("  Step 1: 处理 HRS (美国)\n")
cat("================================================================\n")
hrs_raw <- read_dta(HRS_FILE(),
  col_select = c("hhidpn", paste0("h", 1:15, "atotb"), paste0("h", 1:15, "itot"), "rabplace", "rafeduc", "rameduc"))

hrs_long <- hrs_raw %>%
  pivot_longer(cols=matches("^h[0-9]+(atotb|itot)$"), 
               names_to=c("wave",".value"), names_pattern="h(\\d+)(atotb|itot)") %>%
  mutate(wave = as.integer(wave), id = paste0("HRS_", hhidpn)) %>%
  filter(!is.na(atotb) | !is.na(itot))

hrs_baseline <- hrs_long %>%
  inner_join(baseline_waves %>% filter(grepl("HRS_", id)), by=c("id", "wave"="base_wave")) %>%
  select(id, wealth_total = atotb, income_total = itot, 
         child_dad_edu = rafeduc, child_mom_edu = rameduc, child_hardship = rabplace) %>%
  mutate(study = "USA", country = "USA")

cat("================================================================\n")
cat("  Step 2: 处理 ELSA (英国)\n")
cat("================================================================\n")
elsa_raw <- read_dta(ELSA_FILE(),
  col_select = c("idauniq", paste0("h", 1:9, "atotb"), paste0("h", 1:9, "itot"), "rasfnhch", "ramomeduage", "radadeduage"))

elsa_long <- elsa_raw %>%
  pivot_longer(cols=matches("^h[0-9]+(atotb|itot)$"), 
               names_to=c("wave",".value"), names_pattern="h(\\d+)(atotb|itot)") %>%
  mutate(wave = as.integer(wave), id = paste0("ELSA_", idauniq)) %>%
  filter(!is.na(atotb) | !is.na(itot))

elsa_baseline <- elsa_long %>%
  inner_join(baseline_waves %>% filter(grepl("ELSA_", id)), by=c("id", "wave"="base_wave")) %>%
  select(id, wealth_total = atotb, income_total = itot, 
         child_dad_edu = radadeduage, child_mom_edu = ramomeduage, child_hardship = rasfnhch) %>%
  mutate(study = "England", country = "England")

cat("================================================================\n")
cat("  Step 3: 处理 SHARE (欧洲)\n")
cat("================================================================\n")
share_raw <- read_dta(SHARE_FILE(),
  col_select = c("mergeid", "country", paste0("h", c(1,2,4,5,6,7,8), "atotb"), "h1itot", paste0("h", c(2,4,5,6,7,8), "ittot"), "radadeducl", "ramomeducl"))

names(share_raw)[names(share_raw) == "h1itot"] <- "h1ittot"

share_long <- share_raw %>%
  pivot_longer(cols=matches("^h[0-9]+(atotb|ittot)$"), 
               names_to=c("wave",".value"), names_pattern="h(\\d+)(atotb|ittot)") %>%
  mutate(wave = as.integer(wave), id = paste0("SHARE_", mergeid), country = as.character(haven::as_factor(country))) %>%
  filter(!is.na(atotb) | !is.na(ittot))

share_baseline <- share_long %>%
  inner_join(baseline_waves %>% filter(grepl("SHARE_", id)), by=c("id", "wave"="base_wave")) %>%
  select(id, wealth_total = atotb, income_total = ittot, 
         child_dad_edu = radadeducl, child_mom_edu = ramomeducl, country) %>%
  mutate(study = "Europe", child_hardship = NA)

cat("================================================================\n")
cat("  Step 4: 处理 CHARLS (中国)\n")
cat("================================================================\n")
charls_raw <- read_dta(CHARLS_FILE(),
  col_select = c("ID", "hh1atotb", "hh2atotb", "h3atotb", "h4atotb", 
                 "hh1itot", "hh2itot", "hh3itot", "hh4itot",
                 "radadeducl", "ramomeducl"))

# 修复波次前缀不一致导致的宽转长错位
names(charls_raw)[names(charls_raw) == "h3atotb"] <- "hh3atotb"
names(charls_raw)[names(charls_raw) == "h4atotb"] <- "hh4atotb"

charls_long <- charls_raw %>%
  pivot_longer(cols=matches("^hh[0-9]+(atotb|itot)$"), 
               names_to=c("wave",".value"), names_pattern="hh(\\d+)(atotb|itot)") %>%
  mutate(wave = as.integer(wave), id = paste0("CHARLS_", ID)) %>%
  filter(!is.na(atotb) | !is.na(itot))

charls_baseline <- charls_long %>%
  inner_join(baseline_waves %>% filter(grepl("CHARLS_", id)), by=c("id", "wave"="base_wave")) %>%
  select(id, wealth_total = atotb, income_total = itot, 
         child_dad_edu = radadeducl, child_mom_edu = ramomeducl) %>%
  mutate(study = "China", country = "China", child_hardship = NA)

cat("================================================================\n")
cat("  合并与计算\n")
cat("================================================================\n")

# 移除 labelled 格式以防止 bind_rows 报错
hrs_baseline <- haven::zap_labels(hrs_baseline)
elsa_baseline <- haven::zap_labels(elsa_baseline)
share_baseline <- haven::zap_labels(share_baseline)
charls_baseline <- haven::zap_labels(charls_baseline)

all_ses <- bind_rows(hrs_baseline, elsa_baseline, share_baseline, charls_baseline)

all_ses <- df_main %>% distinct(id, enroll_year, study, female, edu3_f) %>%
  inner_join(all_ses, by=c("id", "study"))

# 过滤 SHARE 2017 入组波次，由于该组财富缺失严重，不参与主干财富分析
all_ses <- all_ses %>%
  mutate(structural_miss = ifelse(study == "Europe" & enroll_year == 2017, 1, 0)) %>%
  filter(structural_miss == 0)

# 国家-波次内部百分位排名
all_ses <- all_ses %>%
  group_by(country, enroll_year) %>%
  mutate(
    wealth_percentile = percent_rank(wealth_total),
    income_percentile = percent_rank(income_total)
  ) %>%
  ungroup()

saveRDS(all_ses, file.path(out_dir, "baseline_ses_wealth.rds"))
cat("财富与SES数据已成功提取，已排除结构性缺失，存入 baseline_ses_wealth.rds\n")
