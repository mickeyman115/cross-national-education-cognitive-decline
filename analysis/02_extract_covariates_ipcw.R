##############################################################################
# 10. 提取协变量与生存状态 (10_extract_covariates_ipw.R)
# 目标：从原始数据库提取婚姻、高血压、糖尿病、心脏病、中风，及访谈状态
##############################################################################

library(haven)
library(dplyr)
library(tidyr)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
out_dir <- PROJECT_DIR
df_clean <- readRDS(file.path(out_dir, "clean_longitudinal_data.rds"))

cat("================================================================\n")
cat("  Step 1: 处理 HRS (美国)\n")
cat("================================================================\n")
hrs_raw <- read_dta(HRS_FILE(),
  col_select = c("hhidpn", 
                 paste0("r",8:13,"iwstat"), paste0("r",8:13,"mstat"),
                 paste0("r",8:13,"diabe"), paste0("r",8:13,"hibpe"), 
                 paste0("r",8:13,"hearte"), paste0("r",8:13,"stroke")))

hrs_long <- hrs_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(iwstat|mstat|diabe|hibpe|hearte|stroke)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(iwstat|mstat|diabe|hibpe|hearte|stroke)") %>%
  mutate(wave=as.integer(wave), id=paste0("HRS_",hhidpn)) %>%
  select(id, wave, iwstat, mstat, diabe, hibpe, hearte, stroke)

cat("================================================================\n")
cat("  Step 2: 处理 ELSA (英国)\n")
cat("================================================================\n")
elsa_raw <- read_dta(ELSA_FILE(),
  col_select = c("idauniq",
                 paste0("r",1:9,"iwstat"), paste0("r",1:9,"mstat"),
                 paste0("r",1:9,"diabe"), paste0("r",1:9,"hibpe"), 
                 paste0("r",1:9,"hearte"), paste0("r",1:9,"stroke")))

elsa_long <- elsa_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(iwstat|mstat|diabe|hibpe|hearte|stroke)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(iwstat|mstat|diabe|hibpe|hearte|stroke)") %>%
  mutate(wave=as.integer(wave), id=paste0("ELSA_",idauniq)) %>%
  select(id, wave, iwstat, mstat, diabe, hibpe, hearte, stroke)

cat("================================================================\n")
cat("  Step 3: 处理 SHARE (欧洲)\n")
cat("================================================================\n")
share_raw <- read_dta(SHARE_FILE(),
  col_select = c("mergeid",
                 paste0("r",c(1,2,4:8),"iwstat"), paste0("r",c(1,2,4:8),"mstat"),
                 paste0("r",c(1,2,4:8),"diabe"), paste0("r",c(1,2,4:8),"hibpe"), 
                 paste0("r",c(1,2,4:8),"hearte"), paste0("r",c(1,2,4:8),"stroke")))

share_long <- share_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(iwstat|mstat|diabe|hibpe|hearte|stroke)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(iwstat|mstat|diabe|hibpe|hearte|stroke)") %>%
  mutate(wave=as.integer(wave), id=paste0("SHARE_",mergeid)) %>%
  select(id, wave, iwstat, mstat, diabe, hibpe, hearte, stroke)

cat("================================================================\n")
cat("  Step 4: 处理 CHARLS (中国)\n")
cat("================================================================\n")
charls_raw <- read_dta(CHARLS_FILE(),
  col_select = c("ID",
                 paste0("r",1:4,"iwstat"), paste0("r",1:4,"mstat"),
                 paste0("r",1:4,"diabe"), paste0("r",1:4,"hibpe"), 
                 paste0("r",1:4,"hearte"), paste0("r",1:4,"stroke")))

charls_long <- charls_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(iwstat|mstat|diabe|hibpe|hearte|stroke)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(iwstat|mstat|diabe|hibpe|hearte|stroke)") %>%
  mutate(wave=as.integer(wave), id=paste0("CHARLS_",ID)) %>%
  select(id, wave, iwstat, mstat, diabe, hibpe, hearte, stroke)

cat("================================================================\n")
cat("  合并所有协变量数据\n")
cat("================================================================\n")

all_cov_long <- bind_rows(hrs_long, elsa_long, share_long, charls_long)

# 编码疾病变量 (0=No, 1=Yes, NA=missing)
# 通常 0=No, 1=Yes, .a/.m/.r = missing (in read_dta these become NA)
all_cov_long <- all_cov_long %>%
  mutate(across(c(diabe, hibpe, hearte, stroke), ~ifelse(.x == 1, 1, ifelse(.x == 0, 0, NA))))

# 编码婚姻状态 (1=Married/Partnered, 0=Unpartnered)
# Harmonized code: 1=Married, 2=Married but separated, 3=Partner, 4=Separated, 5=Divorced, 6=Widowed, 7=Never married, 8=Other
all_cov_long <- all_cov_long %>%
  mutate(married = case_when(
    mstat %in% c(1, 3) ~ 1,
    mstat %in% c(2, 4, 5, 6, 7, 8) ~ 0,
    TRUE ~ NA_real_
  ))

# 编码访谈状态 iwstat
# Harmonized code: 1=Respond, 4=Died, 5/6=Non-response
# 有时 0=Not in wave
cat("iwstat frequency:\n")
print(table(all_cov_long$iwstat, useNA="always"))

cat("================================================================\n")
cat("  与主分析样本合并，提取基线特征\n")
cat("================================================================\n")

# 合并到主样本
df_full <- df_clean %>%
  left_join(all_cov_long, by = c("id", "wave"))

# 提取基线协变量 (第一个有效波次)
df_base_cov <- df_full %>%
  group_by(id) %>%
  arrange(wave) %>%
  slice(1) %>%
  select(id, 
         base_married = married,
         base_diabe = diabe,
         base_hibpe = hibpe,
         base_hearte = hearte,
         base_stroke = stroke) %>%
  ungroup()

df_full <- df_full %>%
  left_join(df_base_cov, by = "id")

cat("基线协变量缺失率 (在纵向样本中):\n")
summary(df_full %>% select(base_married, base_diabe, base_hibpe, base_hearte, base_stroke))

saveRDS(df_full, file.path(out_dir, "data_with_covariates.rds"))
saveRDS(all_cov_long, file.path(out_dir, "all_covariates_long.rds"))

cat("数据提取完成并保存至 data_with_covariates.rds 和 all_covariates_long.rds\n")
