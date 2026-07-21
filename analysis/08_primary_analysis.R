##############################################################################
# 99a_final_primary.R — 全样本 M0 + MICE 童年SES + M1 (m=20)
# 最终分析锁定版本
##############################################################################

library(dplyr)
library(mice)
library(lme4)
library(lmerTest)
library(emmeans)
library(parallel)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
proj_dir <- PROJECT_DIR
output_dir <- file.path(proj_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sink(file.path(output_dir, "final_run.log"), split = TRUE)
cat("================================================================\n")
cat("  99a_final_primary.R — Final Primary Analysis\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))

# ================================================================
# Section 1: Unified contrast extraction function
# ================================================================
extract_midlow <- function(mod, mod_name) {
  cat("  Extracting contrasts for:", mod_name, "\n")

  # 1. Within-study Mid-Low at t=4
  emt <- emtrends(mod, pairwise ~ edu3_f | study, var = "time_in_study",
                  at = list(time_in_study = 4))
  cont_df <- as.data.frame(emt$contrasts)

  # emtrends returns "Low - Mid" (alphabetical); negate for "Mid - Low"
  within <- cont_df %>%
    filter(contrast == "Low - Mid") %>%
    mutate(estimate = -estimate, contrast = "Mid - Low")

  within_res <- data.frame(
    model = mod_name,
    study = as.character(within$study),
    contrast = "Within Mid-Low",
    estimate = within$estimate,
    SE = within$SE,
    stringsAsFactors = FALSE
  )

  # 2. Cross-national interaction: D_China - D_comparator
  emt_full <- emtrends(mod, ~ study * edu3_f, var = "time_in_study",
                       at = list(time_in_study = 4))
  int_comp <- as.data.frame(contrast(emt_full, interaction = c("pairwise", "pairwise")))

  # Filter: China in study_pairwise, Low-Mid in edu3_f_pairwise
  cross_raw <- int_comp %>%
    filter(grepl("China", study_pairwise) & edu3_f_pairwise == "Low - Mid")

  cross_res <- data.frame()
  for (i in 1:nrow(cross_raw)) {
    row <- cross_raw[i, ]
    sp <- as.character(row$study_pairwise)
    parts <- trimws(unlist(strsplit(sp, " - ")))
    comp <- setdiff(parts, "China")
    if (length(comp) == 0) next
    comp <- comp[1]

    # Raw interaction = (StudyA - StudyB)(Low - Mid)
    # If "China - X": (China_Low - China_Mid) - (X_Low - X_Mid) = -(D_China - D_X)
    # If "X - China": (X_Low - X_Mid) - (China_Low - China_Mid) = (D_China - D_X)
    if (parts[1] == "China") {
      est <- -row$estimate
    } else {
      est <- row$estimate
    }

    # ALGEBRAIC CHECK
    d_china <- within_res$estimate[within_res$study == "China"]
    d_comp  <- within_res$estimate[within_res$study == comp]
    err <- abs((d_china - d_comp) - est)
    if (err >= 1e-8) {
      stop(paste("ALGEBRAIC CHECK FAILED for", mod_name, sp,
                 ": manual =", d_china - d_comp, "vs emmeans =", est, "err =", err))
    }

    cross_res <- rbind(cross_res, data.frame(
      model = mod_name,
      study = paste0("China vs ", comp),
      contrast = "Cross-national",
      estimate = est,
      SE = row$SE,
      stringsAsFactors = FALSE
    ))
  }

  rbind(within_res, cross_res)
}

# ================================================================
# Section 2: Load data
# ================================================================
cat("Loading data...\n")
df_ipcw <- readRDS(file.path(proj_dir, "data_with_ipcw.rds"))
parent_ses <- readRDS(file.path(proj_dir, "fullsample_childhood_ses.rds"))

df_model <- df_ipcw %>% filter(!is.na(cogtot))
df_model$edu3_f <- factor(df_model$edu3_f, levels = c("Low", "Mid", "High"))
df_model$study  <- factor(df_model$study, levels = c("China", "USA", "England", "Europe"))

n_ids <- length(unique(df_model$id))
n_obs <- nrow(df_model)
cat("Full sample:", n_ids, "people,", n_obs, "observations\n")
stopifnot(n_ids == 144642)
stopifnot(n_obs == 542426)

# ================================================================
# Section 3: Full-sample M0
# ================================================================
cat("\n================================================================\n")
cat("  Section 3: Full-sample M0\n")
cat("================================================================\n")

base_formula <- cogtot ~ age_base_c + enroll_year_c + retest_flag + female +
  poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
  (1 + time_in_study | id)

m0 <- lmer(base_formula, data = df_model, weights = ipcw, control = ctrl)

m0_contrasts <- extract_midlow(m0, "M0")

m0_diag <- data.frame(
  model = "M0", imputation = 0,
  nobs = nobs(m0),
  n_ids = length(unique(m0@frame$id)),
  is_singular = isSingular(m0),
  converged = length(m0@optinfo$conv$lme4$messages) == 0,
  conv_msg = paste(m0@optinfo$conv$lme4$messages, collapse = "; "),
  grad_norm = if (!is.null(m0@optinfo$derivs$gradient)) max(abs(m0@optinfo$derivs$gradient)) else NA_real_
)
cat("M0 nobs:", nobs(m0), "n_ids:", m0_diag$n_ids,
    "converged:", m0_diag$converged, "singular:", m0_diag$is_singular, "\n")

saveRDS(m0, file.path(output_dir, "fullsample_m0_model.rds"))

cat("\nM0 contrasts:\n")
print(m0_contrasts)

# ================================================================
# Section 4: Full-sample MICE for child_ses (fix 89-person issue)
# ================================================================
cat("\n================================================================\n")
cat("  Section 4: MICE child_ses (m=20, only complete predictors)\n")
cat("================================================================\n")

# Build baseline table — ONLY use predictors guaranteed complete
id_aux <- df_model %>%
  group_by(id) %>%
  summarize(
    study = study[1],
    female = female[1],
    edu3_f = as.character(edu3_f[1]),
    base_age = age_base_c[1],
    enroll_year = enroll_year_c[1],
    n_obs = n(),
    max_time = max(time_in_study, na.rm = TRUE),
    base_cogtot = cogtot[which.min(time_in_study)],
    last_cogtot = cogtot[which.max(time_in_study)],
    cogtot_change = ifelse(n() >= 2 & max(time_in_study) > min(time_in_study),
      (cogtot[which.max(time_in_study)] - cogtot[which.min(time_in_study)]) /
        (time_in_study[which.max(time_in_study)] - time_in_study[which.min(time_in_study)]),
      0),  # single-obs or same-time: zero change
    .groups = "drop"
  )

# Get country from baseline_ses_wealth.rds (only source of SHARE country info)
ses_wealth <- readRDS(file.path(proj_dir, "baseline_ses_wealth.rds"))
id_country <- ses_wealth %>%
  select(id, country) %>%
  mutate(country = as.character(country))

# For non-Europe and Europe people not in baseline_ses_wealth, set country = study
id_aux <- id_aux %>%
  left_join(id_country, by = "id") %>%
  mutate(country = case_when(
    study != "Europe" ~ as.character(study),
    !is.na(country) ~ country,
    TRUE ~ "Europe_Unknown"  # SHARE 2017 enrollees without wealth data
  ))

# Merge parent education and encode child_ses
wide_data <- id_aux %>%
  left_join(parent_ses %>% select(id, child_dad_edu, child_mom_edu), by = "id") %>%
  mutate(
    child_ses_raw = ifelse(
      is.na(child_dad_edu) & is.na(child_mom_edu),
      NA_real_,
      pmax(child_dad_edu, child_mom_edu, na.rm = TRUE)
    ),
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
  mutate(child_ses = factor(child_ses, levels = c("1_Low", "2_Mid", "3_High"), ordered = TRUE)) %>%
  select(-child_dad_edu, -child_mom_edu, -child_ses_raw)

cat("\nChild SES distribution by study (pre-MICE):\n")
for (st in c("China", "USA", "England", "Europe")) {
  cat(st, ":\n")
  print(table(wide_data$child_ses[wide_data$study == st], useNA = "always"))
}

# MICE by study
m_imp <- 20L
base_seed <- 20260720L
all_imputed <- list()
all_logged <- list()

for (st_idx in 1:4) {
  st <- c("China", "USA", "England", "Europe")[st_idx]
  cat("\n--- Imputing:", st, "---\n")

  d_st <- wide_data %>% filter(study == st)
  n_miss <- sum(is.na(d_st$child_ses))
  cat("  Need imputation:", n_miss, "/", nrow(d_st),
      "(", round(100 * n_miss / nrow(d_st), 1), "%)\n")

  if (n_miss == 0) {
    for (i in 1:m_imp) {
      all_imputed[[paste0(st, "_", i)]] <- d_st %>%
        select(id, child_ses) %>% mutate(.imp = i)
    }
    next
  }

  # Select approved predictors only
  pred_vars <- c("female", "edu3_f", "base_age", "enroll_year",
                 "n_obs", "max_time", "base_cogtot", "last_cogtot", "cogtot_change")
  if (st == "Europe") pred_vars <- c(pred_vars, "country")

  d_mice <- d_st %>% select(child_ses, all_of(pred_vars))

  # Hard assertion: all predictors complete for people needing imputation
  need_imp <- is.na(d_mice$child_ses)
  for (v in pred_vars) {
    n_na <- sum(is.na(d_mice[[v]][need_imp]))
    if (n_na > 0) {
      stop(paste("ASSERTION FAILED:", v, "has", n_na,
                 "NAs for people needing child_ses imputation in", st))
    }
  }
  cat("  ✓ All predictors complete for imputation targets\n")

  # Factor conversions
  d_mice$edu3_f <- factor(d_mice$edu3_f)
  if ("country" %in% names(d_mice)) d_mice$country <- factor(d_mice$country)

  # Setup MICE
  ini <- mice(d_mice, m = 1, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  meth[] <- ""
  meth["child_ses"] <- "polr"
  pred[] <- 0
  pred["child_ses", ] <- 1
  pred["child_ses", "child_ses"] <- 0

  # Run
  set.seed(base_seed + st_idx)
  imp <- mice(d_mice, m = m_imp, maxit = 20, method = meth,
              predictorMatrix = pred, printFlag = FALSE)

  # Save full mids object
  saveRDS(imp, file.path(output_dir, paste0("mids_", st, ".rds")))

  # LoggedEvents
  n_logged <- if (is.null(imp$loggedEvents)) 0L else nrow(imp$loggedEvents)
  cat("  loggedEvents:", n_logged, "\n")
  if (n_logged > 0) {
    cat("  Event summary:\n")
    print(table(imp$loggedEvents$meth, imp$loggedEvents$dep))
    all_logged[[st]] <- imp$loggedEvents
  }

  # Save chainMean/chainVar and trace plots
  saveRDS(list(chainMean = imp$chainMean, chainVar = imp$chainVar),
          file.path(output_dir, paste0("mice_chain_", st, ".rds")))
  png(file.path(output_dir, paste0("mice_trace_", st, ".png")), width = 800, height = 600)
  tryCatch(plot(imp), error = function(e) cat("  Trace plot error:", e$message, "\n"))
  dev.off()

  # Extract completed data
  for (i in 1:m_imp) {
    comp <- complete(imp, i)
    comp$id <- d_st$id
    all_imputed[[paste0(st, "_", i)]] <- comp %>%
      select(id, child_ses) %>% mutate(.imp = i)
  }
}

# Combine
child_ses_imputed <- bind_rows(all_imputed)
cat("\nTotal imputed rows:", nrow(child_ses_imputed), "\n")
cat("Per imputation:", nrow(child_ses_imputed) / m_imp, "\n")

# HARD ASSERTION: zero residual NAs
for (i in 1:m_imp) {
  n_na <- sum(is.na(child_ses_imputed$child_ses[child_ses_imputed$.imp == i]))
  if (n_na > 0) stop(paste("RESIDUAL NA FOUND: imputation", i, "has", n_na, "NAs"))
}
cat("✓ Zero residual NAs across all", m_imp, "imputations\n")

# Save
saveRDS(child_ses_imputed, file.path(proj_dir, "fullsample_child_ses_m20.rds"))

# ================================================================
# Section 5: Full-sample M1 × 20
# ================================================================
cat("\n================================================================\n")
cat("  Section 5: Full-sample M1 × 20\n")
cat("================================================================\n")

m1_formula <- update(base_formula,
  . ~ . + poly(time_in_study, 2, raw = TRUE) * child_ses * study)

fit_m1 <- function(i) {
  cat("  Running M1 imputation", i, "...\n")

  imp_base <- child_ses_imputed %>% filter(.imp == i) %>% select(id, child_ses)
  d_imp <- df_model %>% inner_join(imp_base, by = "id")

  # Hard assertions
  actual_nrow <- nrow(d_imp)
  actual_nids <- length(unique(d_imp$id))
  actual_na <- sum(is.na(d_imp$child_ses))
  if (actual_nrow != 542426) stop(paste("M1 imp", i, ": nrow =", actual_nrow, "!= 542426"))
  if (actual_nids != 144642) stop(paste("M1 imp", i, ": n_ids =", actual_nids, "!= 144642"))
  if (actual_na > 0) stop(paste("M1 imp", i, ": child_ses has", actual_na, "NAs"))

  mod <- lmer(m1_formula, data = d_imp, weights = ipcw, control = ctrl)

  assign("d_imp", d_imp, envir = .GlobalEnv)
  contr <- extract_midlow(mod, "M1")
  contr$.imp <- i

  diag <- data.frame(
    model = "M1", imputation = i,
    nobs = nobs(mod),
    n_ids = length(unique(mod@frame$id)),
    is_singular = isSingular(mod),
    converged = length(mod@optinfo$conv$lme4$messages) == 0,
    conv_msg = paste(mod@optinfo$conv$lme4$messages, collapse = "; "),
    grad_norm = if (!is.null(mod@optinfo$derivs$gradient)) max(abs(mod@optinfo$derivs$gradient)) else NA_real_
  )

  list(contrasts = contr, diagnostics = diag)
}

m1_results <- lapply(1:m_imp, fit_m1)

m1_all_contr <- bind_rows(lapply(m1_results, `[[`, "contrasts"))
m1_all_diag  <- bind_rows(lapply(m1_results, `[[`, "diagnostics"))

cat("\n--- M1 Diagnostics ---\n")
print(m1_all_diag)

# ================================================================
# Section 6: Rubin pooling for M1
# ================================================================
cat("\n================================================================\n")
cat("  Section 6: Rubin pooling\n")
cat("================================================================\n")

m1_pooled <- m1_all_contr %>%
  group_by(model, study, contrast) %>%
  summarize(
    Q_bar = mean(estimate),
    U_bar = mean(SE^2),
    B = var(estimate),
    .groups = "drop"
  ) %>%
  mutate(
    T_var = U_bar + (1 + 1/m_imp) * B,
    SE = sqrt(T_var),
    df_rubin = (m_imp - 1) * (1 + U_bar / ((1 + 1/m_imp) * B))^2,
    CI_lo = Q_bar - qt(0.975, df_rubin) * SE,
    CI_hi = Q_bar + qt(0.975, df_rubin) * SE,
    p_value = 2 * pt(abs(Q_bar / SE), df_rubin, lower.tail = FALSE),
    FMI = ((1 + 1/m_imp) * B) / T_var,
    MC_error = sqrt(B / m_imp)
  ) %>%
  rename(estimate = Q_bar) %>%
  select(model, study, contrast, estimate, SE, CI_lo, CI_hi, p_value, df_rubin, FMI, MC_error)

# ================================================================
# Section 7: Combine M0 + M1, Holm correction, final output
# ================================================================
cat("\n================================================================\n")
cat("  Section 7: Assembly\n")
cat("================================================================\n")

# M0 results (single fit, use normal approximation for df since no pooling)
m0_final <- m0_contrasts %>%
  mutate(
    CI_lo = estimate - 1.96 * SE,
    CI_hi = estimate + 1.96 * SE,
    p_value = 2 * pnorm(abs(estimate / SE), lower.tail = FALSE),
    df_rubin = NA_real_,
    FMI = NA_real_,
    MC_error = NA_real_
  )

# Combine
primary_results <- bind_rows(m0_final, m1_pooled)

# Holm correction for cross-national comparisons (3 per model)
primary_results$p_holm <- NA_real_
for (mod_name in c("M0", "M1")) {
  cross_idx <- which(primary_results$model == mod_name &
                     primary_results$contrast == "Cross-national")
  if (length(cross_idx) > 0) {
    primary_results$p_holm[cross_idx] <- p.adjust(
      primary_results$p_value[cross_idx], method = "holm")
  }
}

# Add sample info
primary_results$nobs <- nobs(m0)
primary_results$n_ids <- length(unique(m0@frame$id))

cat("\n--- Final Primary Results ---\n")
print(primary_results %>%
  mutate(across(c(estimate, SE, CI_lo, CI_hi), ~ round(., 4)),
         p_value = signif(p_value, 3),
         p_holm = signif(p_holm, 3),
         FMI = round(FMI, 4),
         MC_error = round(MC_error, 5)))

# Save
write.csv(primary_results, file.path(output_dir, "final_primary_results.csv"), row.names = FALSE)
saveRDS(primary_results, file.path(output_dir, "final_primary_results.rds"))

# Diagnostics
all_diags <- rbind(m0_diag, m1_all_diag)
write.csv(all_diags, file.path(output_dir, "final_model_diagnostics.csv"), row.names = FALSE)

# Sample flow
flow <- data.frame(
  step = c("Full sample people", "Full sample observations",
           "M0 nobs", "M0 n_ids",
           "M1 nobs (per imp)", "M1 n_ids (per imp)",
           "MICE m", "MICE maxit", "MICE seed"),
  value = c(144642, 542426,
            nobs(m0), length(unique(m0@frame$id)),
            m1_all_diag$nobs[1], m1_all_diag$n_ids[1],
            m_imp, 20, base_seed)
)
write.csv(flow, file.path(output_dir, "final_sample_flow.csv"), row.names = FALSE)

# MICE logged events summary
if (length(all_logged) > 0) {
  logged_summary <- bind_rows(lapply(names(all_logged), function(st) {
    all_logged[[st]] %>% mutate(study = st)
  }))
  write.csv(logged_summary, file.path(output_dir, "final_mice_loggedEvents.csv"), row.names = FALSE)
}

# SessionInfo
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "final_sessionInfo.txt"))

cat("\n================================================================\n")
cat("  99a complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")
sink()
