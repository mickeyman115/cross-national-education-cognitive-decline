##############################################################################
# 15. 多重插补汇总模型 (15_run_pooled_models.R)
# 正式版 v4：
#   - M0: 基准 (edu × time × study)
#   - M1: + child_ses × time × study
#   - M2: M1 + wealth × time × study
#   - M3: M1 + income × time × study
#   - M4: M1 + wealth × time × study + income × time × study (联合)
#   - 报告总衰减 (vs M0) 和增量衰减 (vs M1)
#   - 检查 M4 共线性 (VIF)
#   - 验证：无 NA、样本一致、收敛、非奇异
##############################################################################

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(parallel)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
out_dir <- PROJECT_DIR

df_ipcw <- readRDS(file.path(out_dir, "data_with_ipcw.rds"))
imputed_wide <- readRDS(file.path(out_dir, "imputed_baseline_wide.rds"))

m_imputations <- max(imputed_wide$.imp)
cat("插补套数:", m_imputations, "\n")

# 分析样本（排除不在插补集中的人）
imp_ids <- unique(imputed_wide$id[imputed_wide$.imp == 1])
df_model <- df_ipcw %>%
  filter(!is.na(cogtot)) %>%
  filter(id %in% imp_ids)

cat("分析样本:", nrow(df_model), "条观测,", length(unique(df_model$id)), "人\n")

cat("================================================================\n")
cat("  Step 1: 在每套插补数据上拟合 M0–M4 并提取效应\n")
cat("================================================================\n")

# emtrends 提取函数
get_contrasts <- function(mod, mod_name) {
  emt <- emtrends(mod, pairwise ~ edu3_f | study, var = "time_in_study",
                  at = list(time_in_study = 4))

  cont_df <- as.data.frame(emt$contrasts)
  mid_low <- cont_df %>%
    filter(contrast == "Low - Mid") %>%
    mutate(estimate = -estimate, contrast = "Mid - Low")

  res <- mid_low %>%
    mutate(model = mod_name) %>%
    select(study, contrast, estimate, SE, df, model)

  # 跨国交互对比
  emt_int <- emtrends(mod, ~ study * edu3_f, var = "time_in_study",
                      at = list(time_in_study = 4))
  int_comp <- contrast(emt_int, interaction = c("pairwise", "pairwise"))
  int_df <- as.data.frame(int_comp)

  int_target <- int_df %>%
    filter(grepl("Low - Mid", edu3_f_pairwise) & grepl("China", study_pairwise))

  int_res <- int_target %>%
    mutate(model = mod_name, study = study_pairwise,
           contrast = "Mid - Low", estimate = -estimate) %>%
    select(study, contrast, estimate, SE, df, model)

  bind_rows(res, int_res)
}

# 提取 SES × time 交互系数
get_ses_time_coefs <- function(mod, mod_name) {
  fe <- fixef(mod)
  fe_names <- names(fe)
  time_ses <- fe_names[grepl("(wealth_percentile|income_percentile)", fe_names) &
                       grepl("time_in_study", fe_names)]
  if (length(time_ses) > 0) {
    data.frame(model = mod_name, term = time_ses, estimate = fe[time_ses],
               row.names = NULL)
  } else {
    data.frame(model = mod_name, term = "NONE", estimate = NA_real_)
  }
}

results_list <- mclapply(1:m_imputations, function(imp) {
  cat("Running imputation", imp, "...\n")

  imp_base <- imputed_wide %>%
    filter(.imp == imp) %>%
    select(id, child_ses, wealth_percentile, income_percentile)

  d_imp <- df_model %>%
    inner_join(imp_base, by = "id")

  # 验证
  n_na_child <- sum(is.na(d_imp$child_ses))
  n_na_wealth <- sum(is.na(d_imp$wealth_percentile))
  n_na_income <- sum(is.na(d_imp$income_percentile))
  sample_n <- nrow(d_imp)

  if (n_na_child > 0 || n_na_wealth > 0 || n_na_income > 0) {
    cat("  *** imp", imp, "残留 NA: child_ses=", n_na_child,
        "wealth=", n_na_wealth, "income=", n_na_income, "***\n")
  }

  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))

  # === M0: Base ===
  m0 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
               poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
               (1 + time_in_study | id),
             data = d_imp, weights = ipcw, control = ctrl)

  # === M1: + Child SES × time × study ===
  m1 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
               poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
               poly(time_in_study, 2, raw = TRUE) * child_ses * study +
               (1 + time_in_study | id),
             data = d_imp, weights = ipcw, control = ctrl)

  # === M2: M1 + Wealth × time × study ===
  m2 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
               poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
               poly(time_in_study, 2, raw = TRUE) * child_ses * study +
               poly(time_in_study, 2, raw = TRUE) * wealth_percentile * study +
               (1 + time_in_study | id),
             data = d_imp, weights = ipcw, control = ctrl)

  # === M3: M1 + Income × time × study ===
  m3 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
               poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
               poly(time_in_study, 2, raw = TRUE) * child_ses * study +
               poly(time_in_study, 2, raw = TRUE) * income_percentile * study +
               (1 + time_in_study | id),
             data = d_imp, weights = ipcw, control = ctrl)

  # === M4: M1 + Wealth + Income jointly × time × study ===
  m4 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
               poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
               poly(time_in_study, 2, raw = TRUE) * child_ses * study +
               poly(time_in_study, 2, raw = TRUE) * wealth_percentile * study +
               poly(time_in_study, 2, raw = TRUE) * income_percentile * study +
               (1 + time_in_study | id),
             data = d_imp, weights = ipcw, control = ctrl)

  # 提取对比
  res_m0 <- get_contrasts(m0, "M0_Base")
  res_m1 <- get_contrasts(m1, "M1_ChildSES")
  res_m2 <- get_contrasts(m2, "M2_Wealth")
  res_m3 <- get_contrasts(m3, "M3_Income")
  res_m4 <- get_contrasts(m4, "M4_Joint")

  res_all <- bind_rows(res_m0, res_m1, res_m2, res_m3, res_m4)
  res_all$.imp <- imp

  # SES × time 交互系数
  coefs_m2 <- get_ses_time_coefs(m2, "M2_Wealth")
  coefs_m3 <- get_ses_time_coefs(m3, "M3_Income")
  coefs_m4 <- get_ses_time_coefs(m4, "M4_Joint")
  ses_coefs <- bind_rows(coefs_m2, coefs_m3, coefs_m4)
  ses_coefs$.imp <- imp

  # 检查 M4 共线性：wealth 和 income 的相关性
  r_wi <- cor(d_imp$wealth_percentile, d_imp$income_percentile, use = "complete.obs")

  # 诊断
  diag_info <- data.frame(
    .imp = imp,
    sample_n = sample_n,
    na_child = n_na_child,
    na_wealth = n_na_wealth,
    na_income = n_na_income,
    m0_conv = length(m0@optinfo$conv$lme4$messages) == 0,
    m1_conv = length(m1@optinfo$conv$lme4$messages) == 0,
    m2_conv = length(m2@optinfo$conv$lme4$messages) == 0,
    m3_conv = length(m3@optinfo$conv$lme4$messages) == 0,
    m4_conv = length(m4@optinfo$conv$lme4$messages) == 0,
    m0_singular = isSingular(m0),
    m1_singular = isSingular(m1),
    m2_singular = isSingular(m2),
    m3_singular = isSingular(m3),
    m4_singular = isSingular(m4),
    wealth_income_cor = r_wi
  )

  list(contrasts = res_all, diagnostics = diag_info, ses_coefs = ses_coefs)
}, mc.cores = 2)

cat("================================================================\n")
cat("  Step 2: 验证诊断\n")
cat("================================================================\n")

all_contrasts <- bind_rows(lapply(results_list, `[[`, "contrasts"))
all_diags <- bind_rows(lapply(results_list, `[[`, "diagnostics"))
all_ses_coefs <- bind_rows(lapply(results_list, `[[`, "ses_coefs"))

cat("\n--- 收敛与奇异性诊断 ---\n")
print(all_diags %>% select(.imp, sample_n, m0_conv:m4_singular))

cat("\n--- 样本量一致性检查 ---\n")
if (length(unique(all_diags$sample_n)) == 1) {
  cat("  ✓ 所有插补的样本量一致:", unique(all_diags$sample_n), "\n")
} else {
  cat("  ✗ 样本量不一致!\n")
  print(table(all_diags$sample_n))
}

cat("\n--- 残留 NA 检查 ---\n")
if (all(all_diags$na_child == 0 & all_diags$na_wealth == 0 & all_diags$na_income == 0)) {
  cat("  ✓ 所有插补数据中无残留 NA\n")
} else {
  cat("  ✗ 仍存在残留 NA!\n")
}

cat("\n--- M4 共线性检查：wealth-income 相关系数 ---\n")
cat("  Mean r:", round(mean(all_diags$wealth_income_cor), 3),
    "  Range:", round(range(all_diags$wealth_income_cor), 3), "\n")

cat("\n--- 财富/收入 × 时间交互系数 ---\n")
ses_coef_summary <- all_ses_coefs %>%
  group_by(model, term) %>%
  summarize(mean_est = mean(estimate, na.rm = TRUE),
            sd_est = sd(estimate, na.rm = TRUE),
            .groups = "drop")
print(ses_coef_summary, n = 50)

cat("================================================================\n")
cat("  Step 3: Rubin 法则合并\n")
cat("================================================================\n")

pooled_results <- all_contrasts %>%
  group_by(model, study, contrast) %>%
  summarize(
    Q_bar = mean(estimate),
    U_bar = mean(SE^2),
    B = var(estimate),
    T_var = U_bar + (1 + 1/m_imputations) * B,
    SE_pooled = sqrt(T_var),
    CI_lo = Q_bar - qt(0.975, df = (m_imputations - 1) * (1 + (U_bar / ((1 + 1/m_imputations) * B)))^2) * SE_pooled,
    CI_hi = Q_bar + qt(0.975, df = (m_imputations - 1) * (1 + (U_bar / ((1 + 1/m_imputations) * B)))^2) * SE_pooled,
    df_pooled = (m_imputations - 1) * (1 + (U_bar / ((1 + 1/m_imputations) * B)))^2,
    t_stat = Q_bar / SE_pooled,
    p_value = 2 * pt(abs(t_stat), df = df_pooled, lower.tail = FALSE),
    FMI = ((1 + 1/m_imputations) * B) / T_var,
    # Monte Carlo error (标准差除以 sqrt(m))
    MC_error = sqrt(B) / sqrt(m_imputations),
    .groups = "drop"
  )

saveRDS(list(pooled = pooled_results, diagnostics = all_diags,
             ses_coefs = all_ses_coefs, raw_contrasts = all_contrasts),
        file.path(out_dir, "rubin_pooled_results.rds"))

cat("\n--- 汇总结果 (含 95% CI 和 FMI) ---\n")
print(pooled_results %>%
        select(model, study, Q_bar, CI_lo, CI_hi, p_value, FMI, MC_error) %>%
        mutate(across(c(Q_bar, CI_lo, CI_hi), ~round(., 4)),
               p_value = signif(p_value, 3),
               FMI = round(FMI, 4),
               MC_error = round(MC_error, 5)),
      n = 50)

cat("\n--- 序贯衰减计算 (China Mid-Low, 含总衰减和增量衰减) ---\n")
ch_m0 <- pooled_results %>% filter(study == "China", contrast == "Mid - Low", model == "M0_Base") %>% pull(Q_bar)
ch_m1 <- pooled_results %>% filter(study == "China", contrast == "Mid - Low", model == "M1_ChildSES") %>% pull(Q_bar)
ch_m2 <- pooled_results %>% filter(study == "China", contrast == "Mid - Low", model == "M2_Wealth") %>% pull(Q_bar)
ch_m3 <- pooled_results %>% filter(study == "China", contrast == "Mid - Low", model == "M3_Income") %>% pull(Q_bar)
ch_m4 <- pooled_results %>% filter(study == "China", contrast == "Mid - Low", model == "M4_Joint") %>% pull(Q_bar)

cat("\n  M0 (Base):", round(ch_m0, 4), "\n")
cat("  M1 (+ Child SES):", round(ch_m1, 4),
    "  总衰减:", round(100*(ch_m0-ch_m1)/ch_m0, 1), "%\n")
cat("  M2 (M1 + Wealth):", round(ch_m2, 4),
    "  总衰减:", round(100*(ch_m0-ch_m2)/ch_m0, 1), "%",
    "  增量衰减 (M1→M2):", round(100*(ch_m1-ch_m2)/ch_m1, 1), "%\n")
cat("  M3 (M1 + Income):", round(ch_m3, 4),
    "  总衰减:", round(100*(ch_m0-ch_m3)/ch_m0, 1), "%",
    "  增量衰减 (M1→M3):", round(100*(ch_m1-ch_m3)/ch_m1, 1), "%\n")
cat("  M4 (M1 + Wealth + Income):", round(ch_m4, 4),
    "  总衰减:", round(100*(ch_m0-ch_m4)/ch_m0, 1), "%",
    "  增量衰减 (M1→M4):", round(100*(ch_m1-ch_m4)/ch_m1, 1), "%\n")

cat("\n多重插补合并分析已完成并保存。\n")
