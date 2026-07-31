###############################################################################
# Primary combined-weight model without first-return weighting
#
# Combines the defensible covariate specifications into one parsimonious model:
# quadratic follow-up time, cohort-specific flexible baseline-age/time
# relationships, sex/time and entry-year/time interactions, and differential
# retest effects by education and cohort. Natural-spline follow-up time remains
# a separate functional-form sensitivity analysis in script 21; multiplying
# both age and time spline bases produced an over-parameterised 139-column model
# with poor precision and is intentionally not promoted to the primary model.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(splines)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

source("analysis/00_config.R")
proj_dir <- PROJECT_DIR
out_dir <- OUTPUT_DIR

df <- readRDS(file.path(proj_dir, "data_with_ipcw.rds")) %>%
  filter(!is.na(cogtot))
df$edu3_f <- factor(df$edu3_f, levels = c("Low", "Mid", "High"))
df$study <- factor(df$study, levels = c("China", "USA", "England", "Europe"))

ctrl <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 300000),
  check.conv.grad = .makeCC(action = "warning", tol = 0.002)
)

cat("Fitting the prespecified demographic-plus-retest primary model\n")
mod <- lmer(
  cogtot ~
    poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
    ns(age_base_c, df = 3) * poly(time_in_study, 2, raw = TRUE) * study +
    enroll_year_c * poly(time_in_study, 2, raw = TRUE) * study +
    female * poly(time_in_study, 2, raw = TRUE) * study +
    retest_flag * edu3_f * study +
    (1 + time_in_study | id),
  data = df,
  weights = ipcw,
  control = ctrl
)

within <- as.data.frame(confint(emtrends(
  mod,
  pairwise ~ edu3_f | study,
  var = "time_in_study",
  at = list(time_in_study = 4, age_base_c = 0, enroll_year_c = 0),
  data = df,
  weights = "proportional",
  lmer.df = "asymptotic"
)$contrasts, adjust = "none")) %>%
  filter(contrast == "Low - Mid") %>%
  transmute(
    model = "Revised_demographic_retest_M0",
    estimand = "instantaneous_rate_year_4",
    comparison = as.character(study),
    estimate = -estimate,
    SE,
    CI_lo = -asymp.UCL,
    CI_hi = -asymp.LCL
  )

tr_all <- emtrends(
  mod,
  ~ study * edu3_f,
  var = "time_in_study",
  at = list(time_in_study = 4, age_base_c = 0, enroll_year_c = 0),
  data = df,
  weights = "proportional",
  lmer.df = "asymptotic"
)
cross_raw <- as.data.frame(confint(
  contrast(tr_all, interaction = c("pairwise", "pairwise")),
  adjust = "none"
))
cross <- cross_raw %>%
  filter(
    grepl("China", study_pairwise),
    edu3_f_pairwise == "Low - Mid"
  ) %>%
  transmute(
    model = "Revised_demographic_retest_M0",
    estimand = "cross_cohort_difference_year_4",
    comparison = paste0(study_pairwise, " [Mid-Low]"),
    estimate = -estimate,
    SE,
    CI_lo = -asymp.UCL,
    CI_hi = -asymp.LCL
  )

em <- emmeans(
  mod,
  ~ time_in_study * retest_flag * edu3_f | study,
  at = list(
    time_in_study = c(0, 4),
    retest_flag = c(0, 1),
    age_base_c = 0,
    enroll_year_c = 0
  ),
  data = df,
  weights = "proportional",
  lmer.df = "asymptotic"
)

# For cumulative change, respect the observed assessment state: baseline uses
# retest_flag=0 and year 4 uses retest_flag=1. The contrast is
# [(Mid_t4,retest1 - Mid_t0,retest0) - (Low_t4,retest1 - Low_t0,retest0)].
grid_one_study <- em@grid %>% filter(study == levels(df$study)[1])
cum_coef <- rep(0, nrow(grid_one_study))
cum_coef[with(grid_one_study, time_in_study == 4 & retest_flag == 1 & edu3_f == "Mid")] <- 1
cum_coef[with(grid_one_study, time_in_study == 0 & retest_flag == 0 & edu3_f == "Mid")] <- -1
cum_coef[with(grid_one_study, time_in_study == 4 & retest_flag == 1 & edu3_f == "Low")] <- -1
cum_coef[with(grid_one_study, time_in_study == 0 & retest_flag == 0 & edu3_f == "Low")] <- 1

change <- as.data.frame(confint(
  contrast(em, method = list("Mid-Low cumulative change, 0-4" = cum_coef), by = "study"),
  adjust = "none"
)) %>%
  transmute(
    model = "Revised_demographic_retest_M0",
    estimand = "average_annual_change_difference_0_to_4",
    comparison = as.character(study),
    estimate = estimate / 4,
    SE = SE / 4,
    CI_lo = asymp.LCL / 4,
    CI_hi = asymp.UCL / 4
  )

results <- bind_rows(within, cross, change)
diagnostics <- data.frame(
  model = "Revised_demographic_retest_M0",
  nobs = nobs(mod),
  n_ids = n_distinct(df$id),
  rank = qr(model.matrix(mod))$rank,
  columns = ncol(model.matrix(mod)),
  singular = isSingular(mod, tol = 1e-4),
  converged = is.null(mod@optinfo$conv$lme4$messages),
  message = if (is.null(mod@optinfo$conv$lme4$messages)) "" else
    paste(mod@optinfo$conv$lme4$messages, collapse = " | ")
)

write.csv(results, file.path(out_dir, "revised_primary_combined_results.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "revised_primary_combined_diagnostics.csv"), row.names = FALSE)

cat("\nResults\n")
print(results)
cat("\nDiagnostics\n")
print(diagnostics)
cat("\nDemographic-plus-retest primary model complete\n")
