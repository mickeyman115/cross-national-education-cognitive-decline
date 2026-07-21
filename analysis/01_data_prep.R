##############################################################################
# 01. 数据清理与提取 (01_data_prep.R)
# 修正教育定义 (ISCED 3-level) 和时间轴 (age_base, time_in_study, retest)
# 强制仅保留重复观测 (N >= 2)
##############################################################################

library(haven)
library(dplyr)
library(tidyr)
library(readr)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
out_dir <- PROJECT_DIR

cat("================================================================\n")
cat("  Step 1: 处理 HRS (美国)\n")
cat("================================================================\n")
hrs_raw <- read_dta(HRS_FILE(),
  col_select = c("hhidpn","ragender","raeduc",
                 paste0("r",8:13,"agey_b"), paste0("r",8:13,"imrc"), paste0("r",8:13,"dlrc")))

hrs_long <- hrs_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(agey_b|imrc|dlrc)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(agey_b|imrc|dlrc)") %>%
  mutate(wave=as.integer(wave),
         year=case_when(wave==8~2006,wave==9~2008,wave==10~2010,wave==11~2012,wave==12~2014,wave==13~2016,TRUE~NA_real_),
         cogtot=imrc+dlrc, study="USA", id=paste0("HRS_",hhidpn),
         female=as.numeric(ragender==2)) %>%
  rename(age=agey_b) %>%
  # 映射 raeduc 为 ISCED 3-level
  mutate(edu3 = case_when(
    raeduc == 1 ~ 1,                # lt high-school -> Low
    raeduc %in% c(2,3) ~ 2,         # ged, high-school grad -> Mid
    raeduc %in% c(4,5) ~ 3,         # some college, college+ -> High
    TRUE ~ NA_real_
  )) %>%
  select(id, study, wave, year, age, female, edu3, cogtot)

cat("================================================================\n")
cat("  Step 2: 处理 ELSA (英国)\n")
cat("================================================================\n")
elsa_raw <- read_dta(ELSA_FILE(),
  col_select = c("idauniq","ragender","raeducl",
                 paste0("r",1:9,"agey"), paste0("r",1:9,"imrc"), paste0("r",1:9,"dlrc")))
elsa_long <- elsa_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(agey|imrc|dlrc)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(agey|imrc|dlrc)") %>%
  mutate(wave=as.integer(wave),
         year=case_when(wave==1~2002,wave==2~2004,wave==3~2006,wave==4~2008,wave==5~2010,wave==6~2012,wave==7~2014,wave==8~2016,wave==9~2018,TRUE~NA_real_),
         cogtot=imrc+dlrc, study="England", id=paste0("ELSA_",idauniq),
         female=as.numeric(ragender==2),
         edu3 = as.numeric(raeducl)) %>% # raeducl 已经是 1=low, 2=mid, 3=high
  rename(age=agey) %>%
  select(id, study, wave, year, age, female, edu3, cogtot)

cat("================================================================\n")
cat("  Step 3: 处理 SHARE (欧洲)\n")
cat("================================================================\n")
share_raw <- read_dta(SHARE_FILE(),
  col_select = c("mergeid","country","ragender","raeducl",
                 paste0("r",c(1,2,4:8),"agey"), paste0("r",c(1,2,4:8),"imrc"), paste0("r",c(1,2,4:8),"dlrc")))
share_long <- share_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(agey|imrc|dlrc)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(agey|imrc|dlrc)") %>%
  mutate(wave=as.integer(wave),
         year=case_when(wave==1~2004,wave==2~2006,wave==4~2011,wave==5~2013,wave==6~2015,wave==7~2017,wave==8~2020,TRUE~NA_real_),
         cogtot=imrc+dlrc, study="Europe", id=paste0("SHARE_",mergeid),
         female=as.numeric(ragender==2),
         edu3 = as.numeric(raeducl)) %>% # raeducl 已经是 1,2,3
  rename(age=agey) %>%
  select(id, study, wave, year, age, female, edu3, cogtot)

cat("================================================================\n")
cat("  Step 4: 处理 CHARLS (中国)\n")
cat("================================================================\n")
charls_raw <- read_dta(CHARLS_FILE(),
  col_select = c("ID","ragender","raeducl",
                 paste0("r",1:4,"agey"), paste0("r",1:4,"imrc"), paste0("r",1:4,"dlrc")))
charls_long <- charls_raw %>%
  pivot_longer(cols=matches("^r[0-9]+(agey|imrc|dlrc)$"),
               names_to=c("wave",".value"), names_pattern="r(\\d+)(agey|imrc|dlrc)") %>%
  mutate(wave=as.integer(wave),
         year=case_when(wave==1~2011,wave==2~2013,wave==3~2015,wave==4~2018,TRUE~NA_real_),
         cogtot=imrc+dlrc, study="China", id=paste0("CHARLS_",ID),
         female=as.numeric(ragender==2),
         edu3 = as.numeric(raeducl)) %>% # raeducl 已经是 1,2,3
  rename(age=agey) %>%
  select(id, study, wave, year, age, female, edu3, cogtot)


cat("================================================================\n")
cat("  Step 5: 合并与严格纳排 (N >= 2) + 生成时间维度\n")
cat("================================================================\n")

df <- bind_rows(hrs_long, elsa_long, share_long, charls_long) %>%
  filter(!is.na(age), !is.na(cogtot), age >= 50, age <= 100, !is.na(edu3), edu3 %in% 1:3, !is.na(female))

df <- df %>%
  group_by(id) %>%
  arrange(year, wave, age) %>%
  # 计算每个人有效的观测次数
  mutate(n_obs = n()) %>%
  # 强制保留至少 2 次测量的个体
  filter(n_obs >= 2) %>%
  # 提取基线参数和纵向随访参数
  mutate(
    baseline_wave = min(wave),
    enroll_year = min(year),
    age_base = first(age),
    time_in_study = year - first(year),
    retest_flag = ifelse(wave > baseline_wave, 1, 0)
  ) %>%
  ungroup() %>%
  # 定义类别变量
  mutate(
    edu3_f = factor(edu3, levels=c(1,2,3), labels=c("Low", "Mid", "High")),
    study = factor(study, levels=c("China","USA","England","Europe")),
    # 为了模型截距更有实际意义，对基线年龄进行群体中心化 (e.g., 65)
    age_base_c = age_base - 65,
    enroll_year_c = enroll_year - 2011 # 以CHARLS起点为基准
  )

cat("\n合并后严格纵向样本分布 (N obs >= 2):\n")
table(df$study)

# 保存处理好的数据
saveRDS(df, file.path(out_dir, "clean_longitudinal_data.rds"))
cat("数据已保存至 clean_longitudinal_data.rds\n")
