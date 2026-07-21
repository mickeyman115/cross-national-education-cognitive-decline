##############################################################################
# 14. 财富与童年 SES 的多重插补 (14_impute_mice.R)
# 正式版 v4 (m=20)：
#   - 排除 SHARE 2017 (结构性财富缺失)
#   - inner_join 排除无基线财富数据者
#   - ELSA 父母教育 1-7 类别码正确编码
#   - 修复 pmax(NA,NA,na.rm=TRUE) → -Inf 问题
#   - MICE 方法按变量类型指定 (pmm/polr/logreg)
#   - 固定变量不插补
#   - 保存完整 MICE 诊断（mids 对象、loggedEvents、轨迹图）
#   - 纳入分波次认知辅助信息
##############################################################################

library(dplyr)
library(tidyr)
library(mice)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
out_dir <- PROJECT_DIR

df <- readRDS(file.path(out_dir, "data_with_covariates.rds"))
ses_base <- readRDS(file.path(out_dir, "baseline_ses_wealth.rds"))

cat("================================================================\n")
cat("  Step 1: 提取个体级别的基线辅助变量（含认知轨迹摘要）\n")
cat("================================================================\n")

# 计算个体的辅助信息，包括纵向认知轨迹摘要
id_aux <- df %>%
  filter(!is.na(cogtot)) %>%
  group_by(id) %>%
  summarize(
    max_time_in_study = max(time_in_study, na.rm = TRUE),
    n_obs = n(),
    base_cogtot = cogtot[wave == baseline_wave][1],
    # 纵向认知辅助：最后一次认知、认知变化率
    last_cogtot = cogtot[which.max(time_in_study)],
    cogtot_change = ifelse(n() >= 2,
                           (cogtot[which.max(time_in_study)] - cogtot[which.min(time_in_study)]) /
                             (time_in_study[which.max(time_in_study)] - time_in_study[which.min(time_in_study)] + 0.01),
                           NA_real_),
    base_age = age[wave == baseline_wave][1],
    base_married = base_married[wave == baseline_wave][1],
    base_diabe = base_diabe[wave == baseline_wave][1],
    base_hibpe = base_hibpe[wave == baseline_wave][1],
    base_hearte = base_hearte[wave == baseline_wave][1],
    base_stroke = base_stroke[wave == baseline_wave][1],
    female = female[1],
    edu3_f = edu3_f[1],
    study = study[1],
    enroll_year = enroll_year[1],
    .groups = "drop"
  )

# 排除 SHARE 2017
n_before <- nrow(id_aux)
id_aux <- id_aux %>% filter(!(study == "Europe" & enroll_year == 2017))
n_after <- nrow(id_aux)
cat("排除 SHARE 2017 入组者:", n_before - n_after, "人\n")

# inner_join 排除无基线财富数据者
wide_data <- id_aux %>%
  inner_join(ses_base %>% select(id, child_dad_edu, child_mom_edu, child_hardship,
                                 wealth_percentile, income_percentile, country),
            by = "id")

n_dropped <- nrow(id_aux) - nrow(wide_data)
cat("排除无基线财富数据者:", n_dropped, "人\n")
cat("财富分析最终样本:", nrow(wide_data), "人\n")

cat("================================================================\n")
cat("  Step 2: 协调童年 SES (Childhood SES Harmonization)\n")
cat("================================================================\n")

# 安全计算 child_ses_raw，双亲均 NA 时保留 NA
wide_data <- wide_data %>%
  mutate(
    child_ses_raw = ifelse(
      is.na(child_dad_edu) & is.na(child_mom_edu),
      NA_real_,
      pmax(child_dad_edu, child_mom_edu, na.rm = TRUE)
    )
  )

# ELSA: 1-7 类别码 → 1-2=Low, 3-4=Mid, 5-7=High
# HRS: 年限 → <12=Low, 12=Mid, >12=High
# CHARLS/SHARE: 0-3 类别码 → ≤1=Low, 2=Mid, ≥3=High
wide_data <- wide_data %>%
  mutate(
    child_ses = case_when(
      is.na(child_ses_raw) ~ NA_character_,
      study == "USA" & child_ses_raw < 12 ~ "1_Low",
      study == "USA" & child_ses_raw == 12 ~ "2_Mid",
      study == "USA" & child_ses_raw > 12 ~ "3_High",
      study %in% c("China", "Europe") & child_ses_raw <= 1 ~ "1_Low",
      study %in% c("China", "Europe") & child_ses_raw == 2 ~ "2_Mid",
      study %in% c("China", "Europe") & child_ses_raw >= 3 ~ "3_High",
      study == "England" & child_ses_raw <= 2 ~ "1_Low",
      study == "England" & child_ses_raw %in% c(3, 4) ~ "2_Mid",
      study == "England" & child_ses_raw >= 5 ~ "3_High",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(child_ses = factor(child_ses, levels = c("1_Low", "2_Mid", "3_High"), ordered = TRUE))

cat("\n--- 童年 SES 分布诊断 ---\n")
for (st in c("USA", "England", "China", "Europe")) {
  cat("\n", st, ":\n")
  print(table(wide_data$child_ses[wide_data$study == st], useNA = "always"))
}

# 去除原始变量
wide_data <- wide_data %>%
  select(-child_dad_edu, -child_mom_edu, -child_hardship, -child_ses_raw)

cat("================================================================\n")
cat("  Step 3: 按队列执行 MICE 多重插补 (m=20 正式版)\n")
cat("================================================================\n")

m_imputations <- 20  # 正式分析

set.seed(2026)
imputed_list <- list()
mice_diagnostics <- list()

for (st in unique(wide_data$study)) {
  cat("\nImputing baseline data for:", st, "\n")
  d_st <- wide_data %>% filter(study == st) %>% as.data.frame()

  d_imp <- d_st %>% select(-id, -study)

  # 因子化二分类变量
  binary_vars <- c("base_married", "base_diabe", "base_hibpe", "base_hearte", "base_stroke")
  for (v in binary_vars) {
    if (v %in% names(d_imp) && !all(is.na(d_imp[[v]]))) {
      d_imp[[v]] <- as.factor(d_imp[[v]])
    }
  }

  if ("country" %in% names(d_imp)) d_imp$country <- as.factor(d_imp$country)
  d_imp$edu3_f <- as.factor(d_imp$edu3_f)
  d_imp$female <- as.numeric(d_imp$female)

  # MICE 方法指定
  ini <- mice(d_imp, m = 1, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix

  # 固定变量不插补
  no_impute_vars <- c("female", "edu3_f", "enroll_year", "n_obs", "max_time_in_study",
                       "base_age", "country", "base_cogtot", "last_cogtot", "cogtot_change")
  for (v in no_impute_vars) {
    if (v %in% names(meth)) meth[v] <- ""
  }

  # 连续变量 → pmm
  continuous_vars <- c("wealth_percentile", "income_percentile")
  for (v in continuous_vars) {
    if (v %in% names(meth) && meth[v] != "") meth[v] <- "pmm"
  }

  # 有序因子 → polr
  if ("child_ses" %in% names(meth) && meth["child_ses"] != "") {
    meth["child_ses"] <- "polr"
  }

  # 二分类 → logreg
  for (v in binary_vars) {
    if (v %in% names(meth) && meth[v] != "") meth[v] <- "logreg"
  }

  cat("  插补方法:\n")
  print(meth[meth != ""])

  # 运行 MICE (maxit=20 for proper convergence)
  imp_st <- mice(d_imp, m = m_imputations, maxit = 20, method = meth,
                 predictorMatrix = pred, seed = 2026 + which(unique(wide_data$study) == st),
                 printFlag = FALSE)

  # 保存诊断信息
  mice_diagnostics[[st]] <- list(
    logged_events = imp_st$loggedEvents,
    convergence = imp_st$chainMean,
    method = meth[meth != ""],
    n_iter = imp_st$iteration
  )

  cat("  loggedEvents 数量:", nrow(imp_st$loggedEvents), "\n")

  # 保存轨迹图
  tryCatch({
    png(file.path(out_dir, paste0("mice_trace_", st, ".png")), width=1200, height=800)
    plot(imp_st)
    dev.off()
    cat("  轨迹图已保存\n")
  }, error = function(e) cat("  轨迹图保存失败:", e$message, "\n"))

  # 提取 m 套数据
  st_imp_list <- list()
  for (i in 1:m_imputations) {
    comp <- complete(imp_st, i)
    comp$id <- d_st$id
    comp$study <- st
    comp$.imp <- i
    st_imp_list[[i]] <- comp
  }

  imputed_list[[st]] <- bind_rows(st_imp_list)
}

final_imputed_wide <- bind_rows(imputed_list)

cat("\n================================================================\n")
cat("  Step 4: 插补后验证\n")
cat("================================================================\n")

key_vars <- c("child_ses", "wealth_percentile", "income_percentile")
any_na <- FALSE
for (imp_i in 1:m_imputations) {
  d_check <- final_imputed_wide %>% filter(.imp == imp_i)
  for (v in key_vars) {
    n_na <- sum(is.na(d_check[[v]]))
    if (n_na > 0) {
      cat("  *** 警告: imp =", imp_i, ", 变量", v, "仍有", n_na, "个 NA ***\n")
      any_na <- TRUE
    }
  }
}
if (!any_na) cat("  ✓ 所有插补数据中无残留 NA\n")

cat("  每套插补数据的样本量:", final_imputed_wide %>% filter(.imp == 1) %>% nrow(), "\n")
cat("  总插补套数:", m_imputations, "\n")

# 观察值 vs 插补值分布对比
cat("\n--- 观察值 vs 插补值分布 (child_ses, imp=1 vs original) ---\n")
for (st in unique(final_imputed_wide$study)) {
  cat("\n", st, ":\n")
  orig <- wide_data %>% filter(study == st)
  imp1 <- final_imputed_wide %>% filter(study == st, .imp == 1)
  cat("  原始 (含 NA):\n")
  print(prop.table(table(orig$child_ses, useNA = "always")))
  cat("  插补后:\n")
  print(prop.table(table(imp1$child_ses, useNA = "always")))
}

saveRDS(final_imputed_wide, file.path(out_dir, "imputed_baseline_wide.rds"))
saveRDS(mice_diagnostics, file.path(out_dir, "mice_diagnostics.rds"))

cat("\n================================================================\n")
cat("  基线多重插补 (m =", m_imputations, ") 已完成\n")
cat("  imputed_baseline_wide.rds + mice_diagnostics.rds 已保存\n")
cat("================================================================\n")
