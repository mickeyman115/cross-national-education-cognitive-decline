###############################################################################
# 25. Consolidate revised primary and bounded sensitivity outputs
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(splines)
  library(emmeans)
})

source("analysis/00_config.R")
proj_dir <- PROJECT_DIR
out_dir <- OUTPUT_DIR

mod <- readRDS(file.path(out_dir, "first_return_weighted_model.rds"))
d <- readRDS(file.path(proj_dir, "data_with_ipcw.rds")) %>%
  filter(!is.na(cogtot))
d$edu3_f <- factor(d$edu3_f, levels = c("Low", "Mid", "High"))
d$study <- factor(d$study, levels = c("China", "USA", "England", "Europe"))

em <- emmeans(
  mod, ~ time_in_study * retest_flag * edu3_f | study,
  at = list(time_in_study = c(0, 4), retest_flag = c(0, 1),
            age_base_c = 0, enroll_year_c = 0),
  data = d, weights = "proportional", lmer.df = "asymptotic"
)
g <- em@grid %>% filter(study == levels(d$study)[1])
k <- rep(0, nrow(g))
k[with(g, time_in_study == 4 & retest_flag == 1 & edu3_f == "Mid")] <- 1
k[with(g, time_in_study == 0 & retest_flag == 0 & edu3_f == "Mid")] <- -1
k[with(g, time_in_study == 4 & retest_flag == 1 & edu3_f == "Low")] <- -1
k[with(g, time_in_study == 0 & retest_flag == 0 & edu3_f == "Low")] <- 1

change_em <- contrast(
  em, method = list("Mid-Low cumulative change, 0-4" = k), by = "study"
)
within_avg <- as.data.frame(summary(change_em, infer = c(TRUE, TRUE), adjust = "none")) %>%
  transmute(
    model = "Primary_combined_selection_weighted_M0",
    estimand = "average_annual_change_difference_0_to_4",
    comparison_type = "within_cohort",
    comparison = as.character(study),
    estimate = estimate / 4, SE = SE / 4,
    CI_lo = asymp.LCL / 4, CI_hi = asymp.UCL / 4,
    z = z.ratio, p_value = p.value, p_holm = NA_real_
  )

cross_avg <- as.data.frame(summary(
  pairs(update(change_em, by = NULL), adjust = "none"),
  infer = c(TRUE, TRUE), adjust = "none"
)) %>%
  filter(grepl("China", contrast)) %>%
  transmute(
    model = "Primary_combined_selection_weighted_M0",
    estimand = "cross_cohort_average_annual_difference_0_to_4",
    comparison_type = "cross_cohort",
    comparison = contrast,
    estimate = estimate / 4, SE = SE / 4,
    CI_lo = asymp.LCL / 4, CI_hi = asymp.UCL / 4,
    z = z.ratio, p_value = p.value,
    p_holm = p.adjust(p.value, method = "holm")
  )

instant_all <- read.csv(file.path(out_dir, "first_return_weighted_results.csv")) %>%
  filter(estimand %in% c("instantaneous_rate_year_4", "cross_cohort_difference_year_4")) %>%
  mutate(
    comparison_type = if_else(grepl("cross_cohort", estimand), "cross_cohort", "within_cohort"),
    z = estimate / SE,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  ) %>%
  group_by(estimand) %>%
  mutate(p_holm = if_else(comparison_type == "cross_cohort",
                          p.adjust(p_value, method = "holm"), NA_real_)) %>%
  ungroup() %>%
  select(model, estimand, comparison_type, comparison, estimate, SE,
         CI_lo, CI_hi, z, p_value, p_holm)

final_results <- bind_rows(within_avg, cross_avg, instant_all)
write.csv(final_results, file.path(out_dir, "final_revised_primary_results.csv"), row.names = FALSE)

locked <- read.csv(file.path(out_dir, "revised_primary_combined_results.csv"))
selected <- read.csv(file.path(out_dir, "first_return_weighted_results.csv"))
robust_compare <- locked %>%
  select(estimand, comparison, revised_without_first_return = estimate) %>%
  inner_join(
    selected %>% select(estimand, comparison, combined_weight_estimate = estimate),
    by = c("estimand", "comparison")
  ) %>%
  mutate(
    absolute_change = combined_weight_estimate - revised_without_first_return,
    relative_change_percent = 100 * absolute_change / abs(revised_without_first_return)
  )
write.csv(robust_compare, file.path(out_dir, "first_return_robustness_comparison.csv"), row.names = FALSE)

cat("Final revised primary results\n")
print(as.data.frame(final_results), row.names = FALSE)
cat("\nFirst-return robustness comparison\n")
print(as.data.frame(robust_compare), row.names = FALSE)
