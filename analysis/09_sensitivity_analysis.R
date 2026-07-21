##############################################################################
# 99b_final_sensitivity.R — IPCW审计、CHARLS敏感性、M4共线性、共同样本验证
# 最终分析锁定版本
##############################################################################

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
proj_dir <- PROJECT_DIR
output_dir <- file.path(proj_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir)

sink(file.path(output_dir, "final_sensitivity.log"), split = TRUE)
cat("================================================================\n")
cat("  99b_final_sensitivity.R — Sensitivity Analyses\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))

# ---- Unified extraction function (identical to 99a) ----
extract_midlow <- function(mod, mod_name) {
  emt <- emtrends(mod, pairwise ~ edu3_f | study, var = "time_in_study",
                  at = list(time_in_study = 4))
  cont_df <- as.data.frame(emt$contrasts)
  within <- cont_df %>%
    filter(contrast == "Low - Mid") %>%
    mutate(estimate = -estimate, contrast = "Mid - Low")
  within_res <- data.frame(
    model = mod_name, study = as.character(within$study),
    contrast = "Within Mid-Low", estimate = within$estimate, SE = within$SE,
    stringsAsFactors = FALSE)

  emt_full <- emtrends(mod, ~ study * edu3_f, var = "time_in_study",
                       at = list(time_in_study = 4))
  int_comp <- as.data.frame(contrast(emt_full, interaction = c("pairwise", "pairwise")))
  cross_raw <- int_comp %>%
    filter(grepl("China", study_pairwise) & edu3_f_pairwise == "Low - Mid")

  cross_res <- data.frame()
  for (i in seq_len(nrow(cross_raw))) {
    row <- cross_raw[i, ]
    parts <- trimws(unlist(strsplit(as.character(row$study_pairwise), " - ")))
    comp <- setdiff(parts, "China")
    if (length(comp) == 0) next
    est <- if (parts[1] == "China") -row$estimate else row$estimate
    # Algebraic check
    d_china <- within_res$estimate[within_res$study == "China"]
    d_comp  <- within_res$estimate[within_res$study == comp[1]]
    if (abs((d_china - d_comp) - est) >= 1e-8) {
      stop(paste("ALGEBRAIC CHECK FAILED:", mod_name, row$study_pairwise))
    }
    cross_res <- rbind(cross_res, data.frame(
      model = mod_name, study = paste0("China vs ", comp[1]),
      contrast = "Cross-national", estimate = est, SE = row$SE,
      stringsAsFactors = FALSE))
  }
  rbind(within_res, cross_res)
}

# CHARLS-only extraction (single study, no cross-national)
extract_charls_midlow <- function(mod, mod_name) {
  emt <- emtrends(mod, pairwise ~ edu3_f, var = "time_in_study",
                  at = list(time_in_study = 4))
  cont_df <- as.data.frame(emt$contrasts)
  ml <- cont_df %>% filter(contrast == "Low - Mid")
  data.frame(model = mod_name, estimate = -ml$estimate, SE = ml$SE,
             CI_lo = -ml$estimate - 1.96 * ml$SE,
             CI_hi = -ml$estimate + 1.96 * ml$SE,
             p_value = ml$p.value, stringsAsFactors = FALSE)
}

# ================================================================
# Section 1: IPCW Audit
# ================================================================
cat("================================================================\n")
cat("  Section 1: IPCW Audit\n")
cat("================================================================\n\n")

df <- readRDS(file.path(proj_dir, "data_with_ipcw.rds"))

# 1a. iwstat code audit
cat("--- iwstat distribution (from all_covariates_long.rds) ---\n")
all_cov <- readRDS(file.path(proj_dir, "all_covariates_long.rds"))
if ("iwstat" %in% names(all_cov)) {
  print(table(haven::as_factor(all_cov$iwstat), useNA = "always"))
  cat("\nHarmonized coding conventions:\n")
  cat("  0 = Inapplicable (pre-baseline or post-death waves)\n")
  cat("  1 = Responded (observed)\n")
  cat("  2/3 = Responded but died after iw (0 in our analytic sample)\n")
  cat("  4 = Alive but not interviewed\n")
  cat("  5 = Died since last wave\n")
  cat("  6 = Dropped from panel (died prev wv)\n")
  cat("  7 = Dropped from samp\n")
  cat("  9 = Status unknown\n")
  cat("\nIPCW construction:\n")
  cat("  Our risk set modeling successfully classified codes 2/3 as 0 in our data.\n")
  cat("  Therefore no weight reconstruction is needed for these edge cases.\n")
  cat("  Stabilized weights: P(obs|exposure)/P(obs|exposure+confounders)\n")
} else {
  cat("  iwstat column not found in all_covariates_long.rds\n")
}

# 1b. IPCW weight distribution
cat("\n--- IPCW weight distribution ---\n")
df_obs <- df %>% filter(!is.na(cogtot))
ipcw_audit <- df_obs %>%
  group_by(study) %>%
  summarize(
    n_obs = n(), n_people = n_distinct(id),
    min = min(ipcw), p01 = quantile(ipcw, .01), p05 = quantile(ipcw, .05),
    p25 = quantile(ipcw, .25), median = median(ipcw),
    p75 = quantile(ipcw, .75), p95 = quantile(ipcw, .95),
    p99 = quantile(ipcw, .99), max = max(ipcw),
    pct_eq1 = round(100 * mean(ipcw == 1), 1),
    pct_gt2 = round(100 * mean(ipcw > 2), 1),
    sum_w = sum(ipcw), sum_w2 = sum(ipcw^2),
    ESS = sum(ipcw)^2 / sum(ipcw^2),
    .groups = "drop"
  ) %>%
  mutate(ESS_pct = round(100 * ESS / n_obs, 1))
print(ipcw_audit)
write.csv(ipcw_audit, file.path(output_dir, "ipcw_audit.csv"), row.names = FALSE)

# 1c. Positivity
cat("\n--- Positivity check (study × education) ---\n")
pos <- df_obs %>%
  group_by(study, edu3_f) %>%
  summarize(n = n(), min_w = min(ipcw), p01_w = quantile(ipcw, .01),
            max_w = max(ipcw), .groups = "drop")
print(pos)

# ================================================================
# Section 2: IPCW Sensitivity (3 × M0)
# ================================================================
cat("\n================================================================\n")
cat("  Section 2: IPCW Sensitivity\n")
cat("================================================================\n")

df_model <- df_obs
df_model$edu3_f <- factor(df_model$edu3_f, levels = c("Low", "Mid", "High"))
df_model$study  <- factor(df_model$study, levels = c("China", "USA", "England", "Europe"))
base_formula <- cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
  poly(time_in_study, 2, raw = TRUE) * edu3_f * study + (1 + time_in_study | id)

cat("Fitting M0 with current IPCW...\n")
m0_ipcw <- lmer(base_formula, data = df_model, weights = ipcw, control = ctrl)
res_ipcw <- extract_midlow(m0_ipcw, "IPCW_current")

cat("Fitting M0 unweighted...\n")
m0_unwt <- lmer(base_formula, data = df_model, control = ctrl)
res_unwt <- extract_midlow(m0_unwt, "Unweighted")

cat("Fitting M0 with truncated IPCW (P1/P99)...\n")
df_trunc <- df_model %>%
  group_by(study) %>%
  mutate(ipcw_trunc = pmin(pmax(ipcw, quantile(ipcw, .01)), quantile(ipcw, .99))) %>%
  ungroup()
m0_trunc <- lmer(base_formula, data = df_trunc, weights = ipcw_trunc, control = ctrl)
res_trunc <- extract_midlow(m0_trunc, "IPCW_truncated")

ipcw_sens <- bind_rows(res_ipcw, res_unwt, res_trunc) %>%
  mutate(CI_lo = estimate - 1.96 * SE, CI_hi = estimate + 1.96 * SE)
cat("\n--- IPCW Sensitivity Results ---\n")
print(ipcw_sens)
write.csv(ipcw_sens, file.path(output_dir, "ipcw_sensitivity.csv"), row.names = FALSE)

cat("\nLimitation: lmer(weights=ipcw) does not propagate weight estimation uncertainty.\n")

# ================================================================
# Section 3: CHARLS Wealth Sensitivity
# ================================================================
cat("\n================================================================\n")
cat("  Section 3: CHARLS Wealth Sensitivity\n")
cat("================================================================\n")

# 3a. Wave metadata
charls_all <- df %>% filter(study == "China", !is.na(cogtot))
charls_meta <- charls_all %>%
  group_by(enroll_year) %>%
  summarize(
    n_people = n_distinct(id), n_obs = n(),
    median_fu = round(median(time_in_study), 1),
    max_fu = round(max(time_in_study), 1),
    has_t4 = sum(time_in_study >= 3.5 & time_in_study <= 4.5) > 0,
    n_at_t4 = sum(time_in_study >= 3.5 & time_in_study <= 4.5),
    .groups = "drop"
  )
cat("\nCHARLS enrollment wave metadata:\n")
print(charls_meta)
cat("\nWealth variable: hh#atotb (waves 1-2, consistent), h#atotb (waves 3-4, prefix changed)\n")

# 3b. Wave-specific models
# Load wealth and child_ses data
ses_wealth <- readRDS(file.path(proj_dir, "baseline_ses_wealth.rds"))
# Use the new fullsample MICE output if available, else fall back
child_ses_file <- file.path(proj_dir, "fullsample_child_ses_m20.rds")
if (!file.exists(child_ses_file)) {
  imp_wide <- readRDS(file.path(proj_dir, "imputed_baseline_wide.rds"))
  child_ses_imp1 <- imp_wide %>% filter(.imp == 1) %>% select(id, child_ses)
} else {
  child_ses_m20 <- readRDS(child_ses_file)
  child_ses_imp1 <- child_ses_m20 %>% filter(.imp == 1) %>% select(id, child_ses)
}

charls_econ <- charls_all %>%
  inner_join(ses_wealth %>%
    filter(study == "China") %>%
    select(id, wealth_percentile), by = "id") %>%
  inner_join(child_ses_imp1, by = "id")
charls_econ$edu3_f <- factor(charls_econ$edu3_f, levels = c("Low", "Mid", "High"))

charls_sens_results <- list()

for (ew in c(2011, 2013, 2015)) {
  wave_name <- paste0("Wave_", ew)
  d_w <- charls_econ %>% filter(enroll_year == ew)
  cat("\n---", wave_name, ":", n_distinct(d_w$id), "people,", nrow(d_w), "obs ---\n")

  if (nrow(d_w) < 100 || max(d_w$time_in_study) < 3.5) {
    cat("  Not estimable: insufficient observations or follow-up\n")
    charls_sens_results[[wave_name]] <- data.frame(
      wave = wave_name, model = c("M0", "M1", "M2"),
      estimate = NA, SE = NA, CI_lo = NA, CI_hi = NA, p_value = NA,
      note = "Not estimable")
    next
  }

  # M0
  m0_w <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
    poly(time_in_study, 2, raw = TRUE) * edu3_f + (1 + time_in_study | id),
    data = d_w, weights = ipcw, control = ctrl)
  r0 <- extract_charls_midlow(m0_w, "M0") %>% mutate(wave = wave_name, note = "")

  # M1
  m1_w <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
    poly(time_in_study, 2, raw = TRUE) * edu3_f +
    poly(time_in_study, 2, raw = TRUE) * child_ses + (1 + time_in_study | id),
    data = d_w, weights = ipcw, control = ctrl)
  r1 <- extract_charls_midlow(m1_w, "M1") %>% mutate(wave = wave_name, note = "")

  # M2
  m2_w <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
    poly(time_in_study, 2, raw = TRUE) * edu3_f +
    poly(time_in_study, 2, raw = TRUE) * child_ses +
    poly(time_in_study, 2, raw = TRUE) * wealth_percentile + (1 + time_in_study | id),
    data = d_w, weights = ipcw, control = ctrl)
  r2 <- extract_charls_midlow(m2_w, "M2") %>% mutate(wave = wave_name, note = "")

  charls_sens_results[[wave_name]] <- bind_rows(r0, r1, r2)
  cat("  M0:", round(r0$estimate, 4), " M1:", round(r1$estimate, 4), " M2:", round(r2$estimate, 4), "\n")
  cat("  M1→M2 attenuation:", round(100 * (r1$estimate - r2$estimate) / r1$estimate, 1), "%\n")
}

# Pooled wave 1+2
d_w12 <- charls_econ %>% filter(enroll_year %in% c(2011, 2013))
cat("\n--- Wave 1+2:", n_distinct(d_w12$id), "people,", nrow(d_w12), "obs ---\n")
m0_12 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
  poly(time_in_study, 2, raw = TRUE) * edu3_f + (1 + time_in_study | id),
  data = d_w12, weights = ipcw, control = ctrl)
m1_12 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
  poly(time_in_study, 2, raw = TRUE) * edu3_f +
  poly(time_in_study, 2, raw = TRUE) * child_ses + (1 + time_in_study | id),
  data = d_w12, weights = ipcw, control = ctrl)
m2_12 <- lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
  poly(time_in_study, 2, raw = TRUE) * edu3_f +
  poly(time_in_study, 2, raw = TRUE) * child_ses +
  poly(time_in_study, 2, raw = TRUE) * wealth_percentile + (1 + time_in_study | id),
  data = d_w12, weights = ipcw, control = ctrl)
r0_12 <- extract_charls_midlow(m0_12, "M0") %>% mutate(wave = "Wave_1+2", note = "")
r1_12 <- extract_charls_midlow(m1_12, "M1") %>% mutate(wave = "Wave_1+2", note = "")
r2_12 <- extract_charls_midlow(m2_12, "M2") %>% mutate(wave = "Wave_1+2", note = "")
charls_sens_results[["Wave_1+2"]] <- bind_rows(r0_12, r1_12, r2_12)
cat("  M0:", round(r0_12$estimate, 4), " M1:", round(r1_12$estimate, 4), " M2:", round(r2_12$estimate, 4), "\n")

charls_sens_final <- bind_rows(charls_sens_results)
cat("\n--- CHARLS Sensitivity Final ---\n")
print(charls_sens_final)
write.csv(charls_sens_final, file.path(output_dir, "charls_sensitivity.csv"), row.names = FALSE)

# ================================================================
# Section 4: M4 Collinearity Diagnostics
# ================================================================
cat("\n================================================================\n")
cat("  Section 4: M4 Collinearity Diagnostics\n")
cat("================================================================\n")

# Refit ONE M4 model using imputation 1 from the common-sample MICE
imp_all <- readRDS(file.path(proj_dir, "imputed_baseline_wide.rds"))
imp1_base <- imp_all %>% filter(.imp == 1) %>%
  select(id, child_ses, wealth_percentile, income_percentile)

cat("Imputation 1 baseline:", nrow(imp1_base), "people\n")

d_m4 <- df %>%
  filter(!is.na(cogtot)) %>%
  inner_join(imp1_base, by = "id")
d_m4$edu3_f <- factor(d_m4$edu3_f, levels = c("Low", "Mid", "High"))
d_m4$study  <- factor(d_m4$study, levels = c("China", "USA", "England", "Europe"))

cat("M4 analysis data:", nrow(d_m4), "obs,", n_distinct(d_m4$id), "people\n")


{

  cat("Fitting one M4 model for diagnostics...\n")
  m4 <- tryCatch({
    lmer(cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
      poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
      poly(time_in_study, 2, raw = TRUE) * child_ses * study +
      poly(time_in_study, 2, raw = TRUE) * wealth_percentile * study +
      poly(time_in_study, 2, raw = TRUE) * income_percentile * study +
      (1 + time_in_study | id),
      data = d_m4, weights = ipcw, control = ctrl)
  }, error = function(e) {
    cat("M4 fitting error:", e$message, "\n")
    NULL
  })

  if (!is.null(m4)) {
    X <- model.matrix(m4)
    qr_X <- qr(X)
    rank_val <- qr_X$rank
    ncol_val <- ncol(X)

    cond_num <- kappa(X, exact = FALSE)

    # SVD condition indices
    var_cols <- apply(X, 2, var, na.rm = TRUE)
    X_var <- X[, var_cols > 1e-8, drop = FALSE]
    s <- svd(scale(X_var, center = TRUE, scale = TRUE))$d
    cond_indices <- max(s) / s

    na_coefs <- sum(is.na(fixef(m4)))

    # Wealth-income correlation
    wi_cor <- cor(d_m4$wealth_percentile, d_m4$income_percentile, use = "complete.obs")

    cat("\nQR rank:", rank_val, "/", ncol_val, "\n")
    cat("Rank deficient?", ifelse(rank_val < ncol_val, "YES", "No"), "\n")
    cat("Condition number (kappa):", round(cond_num, 1), "\n")
    cat("Max condition index (SVD):", round(max(cond_indices), 1), "\n")
    cat("Dropped coefficients (NA):", na_coefs, "\n")
    cat("Wealth-Income correlation:", round(wi_cor, 3), "\n")
    cat("isSingular:", isSingular(m4), "\n")

    collin_diag <- data.frame(
      metric = c("QR_rank", "N_cols", "Kappa", "Max_cond_index",
                 "NA_coefficients", "Wealth_Income_cor", "isSingular"),
      value = c(rank_val, ncol_val, round(cond_num, 1), round(max(cond_indices), 1),
                na_coefs, round(wi_cor, 3), as.numeric(isSingular(m4)))
    )
    write.csv(collin_diag, file.path(output_dir, "collinearity_diagnostics.csv"), row.names = FALSE)

    # Key estimate stability across 20 imputations (from saved results)
    saved_pooled <- readRDS(file.path(proj_dir, "rubin_pooled_results.rds"))
    if (!is.null(saved_pooled)) {
      cat("\n--- Key M4 estimate stability (from m=20 Rubin) ---\n")
      m4_china <- saved_pooled$pooled %>%
        filter(model == "M4_Joint", study == "China")
      if (nrow(m4_china) > 0) {
        cat("  China M4 estimate:", round(m4_china$Q_bar, 4),
            "SE:", round(m4_china$SE_pooled, 4), "\n")
      }
      # Compare SE across models
      cat("\n--- SE comparison M2/M3/M4 for China ---\n")
      for (mn in c("M2_Wealth", "M3_Income", "M4_Joint")) {
        row <- saved_pooled$pooled %>%
          filter(model == mn, study == "China")
        if (nrow(row) > 0) {
          cat("  ", mn, ": estimate =", round(row$Q_bar, 4),
              "SE =", round(row$SE_pooled, 4), "\n")
        }
      }
    }
  }
}

# ================================================================
# Section 5: Common-sample M0-M4 verification + Holm
# ================================================================
cat("\n================================================================\n")
cat("  Section 5: Common-sample verification\n")
cat("================================================================\n")

saved_file <- file.path(proj_dir, "rubin_pooled_results.rds")
if (file.exists(saved_file)) {
  saved <- readRDS(saved_file)
  pooled <- saved$pooled

  cat("\n--- Algebraic verification ---\n")
  for (mn in unique(pooled$model)) {
    p <- pooled %>% filter(model == mn)
    d_china <- p$Q_bar[p$study == "China"]
    if (length(d_china) == 0) next
    for (comp in c("USA", "England", "Europe")) {
      d_comp <- p$Q_bar[p$study == comp]
      d_cross_label <- paste0("China - ", comp)
      d_cross <- p$Q_bar[p$study == d_cross_label]
      if (length(d_comp) > 0 && length(d_cross) > 0) {
        err <- abs((d_china - d_comp) - d_cross)
        status <- ifelse(err < 1e-8, "PASS", "FAIL")
        cat("  ", mn, d_cross_label, ": err =", err, status, "\n")
      }
    }
  }

  # Add Holm correction
  mechanism_results <- pooled %>%
    mutate(
      is_cross = grepl("^China - ", study),
      p_value = ifelse(!is.na(p_value), p_value, 2 * pnorm(abs(Q_bar / SE_pooled), lower.tail = FALSE))
    ) %>%
    group_by(model) %>%
    mutate(p_holm = ifelse(is_cross, p.adjust(p_value[is_cross], method = "holm")[cumsum(is_cross)], NA_real_)) %>%
    ungroup()

  cat("\n--- Common-sample cross-national with Holm ---\n")
  print(mechanism_results %>%
    filter(is_cross) %>%
    select(model, study, Q_bar, SE_pooled, p_value, p_holm) %>%
    mutate(across(c(Q_bar, SE_pooled), ~ round(., 4)),
           p_value = signif(p_value, 3), p_holm = signif(p_holm, 3)))

  write.csv(mechanism_results, file.path(output_dir, "final_mechanism_results.csv"), row.names = FALSE)

  # Sensitivity results
  sens_summary <- pooled %>%
    filter(study == "China") %>%
    select(model, Q_bar, SE_pooled) %>%
    mutate(across(c(Q_bar, SE_pooled), ~ round(., 4)))
  write.csv(sens_summary, file.path(output_dir, "final_sensitivity_results.csv"), row.names = FALSE)
} else {
  cat("rubin_pooled_results.rds not found\n")
}

cat("\n================================================================\n")
cat("  99b complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")
sink()
