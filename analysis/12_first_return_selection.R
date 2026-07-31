###############################################################################
# 23. First-return selection audit and combined-weight sensitivity analysis
#
# Target population: participants eligible at their first valid cognitive
# assessment, including people with only one valid assessment. The outcome
# model remains restricted to repeated observers because a within-person slope
# is unidentified for singletons. A cohort-specific stabilized inverse-
# probability-of-first-return weight reweights repeated observers toward the
# baseline-eligible population; it is multiplied by the existing conditional
# response IPCW. This is a sensitivity analysis under a missing-at-random
# assumption, not proof that attrition bias is absent.
###############################################################################

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(splines)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

source("analysis/00_config.R")
proj_dir <- PROJECT_DIR
out_dir <- OUTPUT_DIR

make_long <- function(path, id_var, prefix, waves, age_suffix, study_name,
                      year_map, education = c("harmonized", "hrs"), include_country = FALSE) {
  education <- match.arg(education)
  id_sym <- rlang::sym(id_var)
  age_vars <- paste0("r", waves, age_suffix)
  imrc_vars <- paste0("r", waves, "imrc")
  dlrc_vars <- paste0("r", waves, "dlrc")
  cols <- c(id_var, "ragender", if (education == "hrs") "raeduc" else "raeducl",
            if (include_country) "country", age_vars, imrc_vars, dlrc_vars)
  x <- read_dta(path, col_select = all_of(cols))
  x %>%
    pivot_longer(
      cols = matches(paste0("^r[0-9]+(", age_suffix, "|imrc|dlrc)$")),
      names_to = c("wave", ".value"),
      names_pattern = paste0("r(\\d+)(", age_suffix, "|imrc|dlrc)")
    ) %>%
    mutate(
      wave = as.integer(wave),
      year = unname(year_map[as.character(wave)]),
      cogtot = imrc + dlrc,
      study = study_name,
      id = paste0(prefix, as.character(!!id_sym)),
      female = as.numeric(ragender == 2),
      edu3 = if (education == "hrs") {
        case_when(
          raeduc == 1 ~ 1,
          raeduc %in% c(2, 3) ~ 2,
          raeduc %in% c(4, 5) ~ 3,
          TRUE ~ NA_real_
        )
      } else as.numeric(raeducl)
    ) %>%
    rename(age = all_of(age_suffix)) %>%
    select(id, study, wave, year, age, female, edu3, cogtot)
}

cat("Reconstructing the baseline-eligible population, including singletons\n")
hrs <- make_long(
  HRS_FILE(),
  "hhidpn", "HRS_", 8:13, "agey_b", "USA",
  c(`8` = 2006, `9` = 2008, `10` = 2010, `11` = 2012, `12` = 2014, `13` = 2016),
  "hrs"
)
elsa <- make_long(
  ELSA_FILE(),
  "idauniq", "ELSA_", 1:9, "agey", "England",
  setNames(seq(2002, 2018, by = 2), 1:9), "harmonized"
)
share <- make_long(
  SHARE_FILE(),
  "mergeid", "SHARE_", c(1, 2, 4:8), "agey", "Europe",
  c(`1` = 2004, `2` = 2006, `4` = 2011, `5` = 2013, `6` = 2015, `7` = 2017, `8` = 2020),
  "harmonized", TRUE
)
charls <- make_long(
  CHARLS_FILE(),
  "ID", "CHARLS_", 1:4, "agey", "China",
  c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018), "harmonized"
)

eligible_long <- bind_rows(hrs, elsa, share, charls) %>%
  filter(
    !is.na(age), !is.na(cogtot), age >= 50, age <= 100,
    !is.na(edu3), edu3 %in% 1:3, !is.na(female)
  ) %>%
  arrange(id, year, wave, age)

eligible_base <- eligible_long %>%
  group_by(id) %>%
  summarise(
    study = first(study),
    baseline_wave = first(wave),
    enroll_year = first(year),
    age_base = first(age),
    cogtot_base = first(cogtot),
    female = first(female),
    edu3 = first(edu3),
    n_valid = n(),
    has_followup = as.integer(n_valid >= 2),
    .groups = "drop"
  ) %>%
  mutate(
    study = factor(study, levels = c("China", "USA", "England", "Europe")),
    edu3_f = factor(edu3, levels = 1:3, labels = c("Low", "Mid", "High")),
    enroll_year_c = enroll_year - 2011
  )

# Baseline health and marital status are taken from the same wave as the first
# valid cognitive assessment. Missing values are retained as explicit factor
# levels rather than deleting baseline-eligible participants.
cov_long <- readRDS(file.path(proj_dir, "all_covariates_long.rds"))
base_cov <- eligible_base %>%
  select(id, baseline_wave) %>%
  left_join(cov_long, by = c("id", "baseline_wave" = "wave")) %>%
  transmute(
    id,
    married = as.character(as.numeric(zap_labels(married))),
    diabe = as.character(as.numeric(zap_labels(diabe))),
    hibpe = as.character(as.numeric(zap_labels(hibpe))),
    hearte = as.character(as.numeric(zap_labels(hearte))),
    stroke = as.character(as.numeric(zap_labels(stroke)))
  ) %>%
  mutate(across(-id, ~factor(replace_na(.x, "Missing"))))

eligible_base <- eligible_base %>% left_join(base_cov, by = "id")

selection_flow <- eligible_base %>%
  group_by(study) %>%
  summarise(
    baseline_eligible = n(),
    repeated_observers = sum(has_followup),
    singletons = sum(1 - has_followup),
    followup_percent = 100 * mean(has_followup),
    .groups = "drop"
  )
write.csv(selection_flow, file.path(out_dir, "first_return_selection_flow.csv"), row.names = FALSE)

selection_profile <- eligible_base %>%
  group_by(study, has_followup) %>%
  summarise(
    n = n(),
    age_mean = mean(age_base),
    female_percent = 100 * mean(female),
    cognition_mean = mean(cogtot_base),
    low_education_percent = 100 * mean(edu3_f == "Low"),
    mid_education_percent = 100 * mean(edu3_f == "Mid"),
    high_education_percent = 100 * mean(edu3_f == "High"),
    .groups = "drop"
  )
write.csv(selection_profile, file.path(out_dir, "first_return_selection_profile.csv"), row.names = FALSE)

cat("Estimating cohort-specific stabilized first-return weights\n")
eligible_base$first_return_sw <- NA_real_
weight_diagnostics <- list()
for (st in levels(eligible_base$study)) {
  d <- eligible_base %>% filter(study == st)
  num <- glm(has_followup ~ edu3_f + female, family = binomial(), data = d)
  den <- glm(
    has_followup ~ edu3_f + female + ns(age_base, df = 3) +
      ns(cogtot_base, df = 3) + enroll_year_c +
      married + diabe + hibpe + hearte + stroke,
    family = binomial(), data = d
  )
  p_num <- pmin(pmax(predict(num, type = "response"), 0.01), 0.99)
  p_den <- pmin(pmax(predict(den, type = "response"), 0.01), 0.99)
  sw <- p_num / p_den
  keep <- d$has_followup == 1
  lo <- quantile(sw[keep], 0.01, na.rm = TRUE)
  hi <- quantile(sw[keep], 0.99, na.rm = TRUE)
  sw <- pmin(pmax(sw, lo), hi)
  eligible_base$first_return_sw[eligible_base$study == st] <- ifelse(keep, sw, NA_real_)
  weight_diagnostics[[st]] <- data.frame(
    study = st,
    numerator_formula = paste(deparse(formula(num)), collapse = ""),
    denominator_formula = paste(deparse(formula(den)), collapse = ""),
    retained_n = sum(keep),
    weight_min = min(sw[keep]),
    weight_p50 = median(sw[keep]),
    weight_p99 = quantile(sw[keep], 0.99),
    weight_max = max(sw[keep]),
    effective_sample_size = sum(sw[keep])^2 / sum(sw[keep]^2),
    c_statistic = suppressWarnings(as.numeric(pROC::auc(d$has_followup, p_den, quiet = TRUE)))
  )
}
weight_diagnostics <- bind_rows(weight_diagnostics)
write.csv(weight_diagnostics, file.path(out_dir, "first_return_weight_diagnostics.csv"), row.names = FALSE)
saveRDS(
  eligible_base %>%
    select(id, study, baseline_wave, enroll_year, has_followup, first_return_sw),
  file.path(out_dir, "first_return_weights_by_id.rds")
)

analysis <- readRDS(file.path(proj_dir, "data_with_ipcw.rds")) %>%
  filter(!is.na(cogtot)) %>%
  left_join(eligible_base %>% select(id, first_return_sw), by = "id") %>%
  mutate(combined_ipcw = ipcw * first_return_sw)
stopifnot(!anyNA(analysis$combined_ipcw))

# Truncate the final product within cohort to retain the cohort-specific weight
# scale and avoid a few compounded extreme weights.
analysis <- analysis %>%
  group_by(study) %>%
  mutate(
    combined_ipcw = pmin(
      pmax(combined_ipcw, quantile(combined_ipcw, 0.01)),
      quantile(combined_ipcw, 0.99)
    )
  ) %>%
  ungroup()
analysis$edu3_f <- factor(analysis$edu3_f, levels = c("Low", "Mid", "High"))
analysis$study <- factor(analysis$study, levels = c("China", "USA", "England", "Europe"))

ctrl <- lmerControl(
  optimizer = "bobyqa", optCtrl = list(maxfun = 300000),
  check.conv.grad = .makeCC(action = "warning", tol = 0.002)
)
cat("Fitting the primary model with combined first-return and conditional-response weights\n")
mod <- lmer(
  cogtot ~
    poly(time_in_study, 2, raw = TRUE) * edu3_f * study +
    ns(age_base_c, df = 3) * poly(time_in_study, 2, raw = TRUE) * study +
    enroll_year_c * poly(time_in_study, 2, raw = TRUE) * study +
    female * poly(time_in_study, 2, raw = TRUE) * study +
    retest_flag * edu3_f * study +
    (1 + time_in_study | id),
  data = analysis, weights = combined_ipcw, control = ctrl
)

within <- as.data.frame(confint(emtrends(
  mod, pairwise ~ edu3_f | study, var = "time_in_study",
  at = list(time_in_study = 4, age_base_c = 0, enroll_year_c = 0),
  data = analysis, weights = "proportional", lmer.df = "asymptotic"
)$contrasts, adjust = "none")) %>%
  filter(contrast == "Low - Mid") %>%
  transmute(
    estimand = "instantaneous_rate_year_4", comparison = as.character(study),
    estimate = -estimate, SE, CI_lo = -asymp.UCL, CI_hi = -asymp.LCL
  )

tr_all <- emtrends(
  mod, ~ study * edu3_f, var = "time_in_study",
  at = list(time_in_study = 4, age_base_c = 0, enroll_year_c = 0),
  data = analysis, weights = "proportional", lmer.df = "asymptotic"
)
cross <- as.data.frame(confint(
  contrast(tr_all, interaction = c("pairwise", "pairwise")), adjust = "none"
)) %>%
  filter(grepl("China", study_pairwise), edu3_f_pairwise == "Low - Mid") %>%
  transmute(
    estimand = "cross_cohort_difference_year_4",
    comparison = paste0(study_pairwise, " [Mid-Low]"),
    estimate = -estimate, SE, CI_lo = -asymp.UCL, CI_hi = -asymp.LCL
  )

em <- emmeans(
  mod, ~ time_in_study * retest_flag * edu3_f | study,
  at = list(time_in_study = c(0, 4), retest_flag = c(0, 1),
            age_base_c = 0, enroll_year_c = 0),
  data = analysis, weights = "proportional", lmer.df = "asymptotic"
)
grid_one <- em@grid %>% filter(study == levels(analysis$study)[1])
k <- rep(0, nrow(grid_one))
k[with(grid_one, time_in_study == 4 & retest_flag == 1 & edu3_f == "Mid")] <- 1
k[with(grid_one, time_in_study == 0 & retest_flag == 0 & edu3_f == "Mid")] <- -1
k[with(grid_one, time_in_study == 4 & retest_flag == 1 & edu3_f == "Low")] <- -1
k[with(grid_one, time_in_study == 0 & retest_flag == 0 & edu3_f == "Low")] <- 1
change <- as.data.frame(confint(
  contrast(em, method = list("Mid-Low cumulative change, 0-4" = k), by = "study"),
  adjust = "none"
)) %>%
  transmute(
    estimand = "average_annual_change_difference_0_to_4",
    comparison = as.character(study), estimate = estimate / 4, SE = SE / 4,
    CI_lo = asymp.LCL / 4, CI_hi = asymp.UCL / 4
  )

results <- bind_rows(within, cross, change) %>%
  mutate(model = "Revised_M0_plus_first_return_and_conditional_IPCW", .before = 1)
diagnostics <- data.frame(
  model = "Revised_M0_plus_first_return_and_conditional_IPCW",
  observations = nobs(mod), people = n_distinct(analysis$id),
  baseline_eligible_people = nrow(eligible_base),
  rank = qr(model.matrix(mod))$rank, columns = ncol(model.matrix(mod)),
  singular = isSingular(mod, tol = 1e-4),
  converged = is.null(mod@optinfo$conv$lme4$messages),
  convergence_message = if (is.null(mod@optinfo$conv$lme4$messages)) "" else
    paste(mod@optinfo$conv$lme4$messages, collapse = " | ")
)

write.csv(results, file.path(out_dir, "first_return_weighted_results.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "first_return_model_diagnostics.csv"), row.names = FALSE)
saveRDS(mod, file.path(out_dir, "first_return_weighted_model.rds"))
cat("\nSelection flow\n"); print(selection_flow)
cat("\nWeight diagnostics\n"); print(weight_diagnostics)
cat("\nResults\n"); print(results)
cat("\nModel diagnostics\n"); print(diagnostics)
cat("\nFirst-return selection analysis complete\n")
