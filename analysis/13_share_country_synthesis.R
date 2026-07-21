###############################################################################
# 24. SHARE country-specific heterogeneity analysis and random-effects synthesis
#
# Every included country uses the same fixed and random-effects specification.
# There is no country-specific fallback model. The analysis is restricted to
# Low and Mid education because that is the primary contrast. Countries require
# >=200 repeated observers, >=50 people in each education group, and >=50
# observations between years 3 and 5 and at least four years of follow-up to
# support the year-4 estimand and identify a random slope.
###############################################################################

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(splines)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
})

source("analysis/00_config.R")
proj_dir <- PROJECT_DIR
out_dir <- OUTPUT_DIR

share_path <- SHARE_FILE()
share_ids <- read_dta(share_path, col_select = c("mergeid", "country"))
country_labels <- attributes(share_ids$country)$labels
country_map <- tibble(
  country_code = as.numeric(country_labels),
  country = names(country_labels)
)
id_country <- share_ids %>%
  transmute(id = paste0("SHARE_", mergeid), country_code = as.numeric(country)) %>%
  left_join(country_map, by = "country_code") %>%
  distinct(id, .keep_all = TRUE)

d <- readRDS(file.path(proj_dir, "data_with_ipcw.rds")) %>%
  filter(study == "Europe", !is.na(cogtot), edu3_f %in% c("Low", "Mid")) %>%
  left_join(
    readRDS(file.path(out_dir, "first_return_weights_by_id.rds")) %>%
      select(id, first_return_sw),
    by = "id"
  ) %>%
  left_join(id_country %>% select(id, country), by = "id") %>%
  filter(!is.na(country)) %>%
  mutate(
    edu2 = factor(edu3_f, levels = c("Low", "Mid")),
    country = as.character(country),
    combined_ipcw = ipcw * first_return_sw
  ) %>%
  mutate(
    combined_ipcw = pmin(
      pmax(combined_ipcw, quantile(combined_ipcw, 0.01)),
      quantile(combined_ipcw, 0.99)
    )
  )

support <- d %>%
  group_by(country) %>%
  summarise(
    people = n_distinct(id), observations = n(),
    low_people = n_distinct(id[edu2 == "Low"]),
    mid_people = n_distinct(id[edu2 == "Mid"]),
    observations_year3_to_5 = sum(time_in_study >= 3 & time_in_study <= 5),
    max_followup = max(time_in_study),
    include = people >= 200 & low_people >= 50 & mid_people >= 50 &
      observations_year3_to_5 >= 50 & max_followup >= 4,
    .groups = "drop"
  ) %>%
  arrange(desc(people))
write.csv(support, file.path(out_dir, "share_country_support.csv"), row.names = FALSE)

ctrl <- lmerControl(
  optimizer = "bobyqa", optCtrl = list(maxfun = 200000),
  check.conv.grad = .makeCC(action = "warning", tol = 0.002)
)
common_formula <- cogtot ~
  poly(time_in_study, 2, raw = TRUE) * edu2 +
  ns(age_base_c, df = 3) * poly(time_in_study, 2, raw = TRUE) +
  female * poly(time_in_study, 2, raw = TRUE) +
  retest_flag * edu2 +
  (1 + time_in_study | id)

extract_country <- function(dd, cname) {
  fit <- tryCatch(
    lmer(common_formula, data = dd, weights = combined_ipcw, control = ctrl),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(result = tibble(country = cname, status = "fit_error", message = fit$message), fit = NULL))
  }
  conv_msg <- fit@optinfo$conv$lme4$messages
  singular <- isSingular(fit, tol = 1e-4)
  status <- if (singular) {
    "singular_random_slope"
  } else if (is.null(conv_msg)) {
    "ok"
  } else {
    "convergence_warning"
  }

  instant <- tryCatch({
    as.data.frame(confint(emtrends(
      fit, pairwise ~ edu2, var = "time_in_study",
      at = list(time_in_study = 4, age_base_c = 0),
      data = dd, weights = "proportional", lmer.df = "asymptotic"
    )$contrasts, adjust = "none")) %>%
      filter(contrast == "Low - Mid") %>%
      transmute(
        country = cname, estimand = "instantaneous_rate_year_4",
        estimate = -estimate, SE, CI_lo = -asymp.UCL, CI_hi = -asymp.LCL
      )
  }, error = function(e) tibble())

  em <- emmeans(
    fit, ~ time_in_study * retest_flag * edu2,
    at = list(time_in_study = c(0, 4), retest_flag = c(0, 1), age_base_c = 0),
    data = dd, weights = "proportional", lmer.df = "asymptotic"
  )
  g <- em@grid
  k <- rep(0, nrow(g))
  k[with(g, time_in_study == 4 & retest_flag == 1 & edu2 == "Mid")] <- 1
  k[with(g, time_in_study == 0 & retest_flag == 0 & edu2 == "Mid")] <- -1
  k[with(g, time_in_study == 4 & retest_flag == 1 & edu2 == "Low")] <- -1
  k[with(g, time_in_study == 0 & retest_flag == 0 & edu2 == "Low")] <- 1
  change <- as.data.frame(confint(
    contrast(em, method = list("Mid-Low cumulative change, 0-4" = k)), adjust = "none"
  )) %>%
    transmute(
      country = cname, estimand = "average_annual_change_difference_0_to_4",
      estimate = estimate / 4, SE = SE / 4,
      CI_lo = asymp.LCL / 4, CI_hi = asymp.UCL / 4
    )

  info <- support %>% filter(country == cname)
  ans <- bind_rows(instant, change) %>%
    mutate(
      people = info$people, observations = info$observations,
      low_people = info$low_people, mid_people = info$mid_people,
      status = status,
      singular = singular,
      rank = qr(model.matrix(fit))$rank,
      columns = ncol(model.matrix(fit)),
      message = if (is.null(conv_msg)) "" else paste(conv_msg, collapse = " | ")
    )
  list(result = ans, fit = fit)
}

included <- support %>% filter(include) %>% pull(country)
cat("Fitting", length(included), "SHARE countries with one common specification\n")
country_results <- list()
for (cname in included) {
  cat("  ", cname, "\n")
  dd <- d %>% filter(country == cname) %>% droplevels()
  z <- extract_country(dd, cname)
  country_results[[cname]] <- z$result
}
country_results <- bind_rows(country_results)
write.csv(country_results, file.path(out_dir, "share_country_common_model_results.csv"), row.names = FALSE)

# Restricted maximum-likelihood random-effects meta-analysis implemented here
# to avoid adding a package dependency. Country estimates are independent
# because each was fitted in a disjoint national sample.
reml_meta <- function(yi, sei) {
  vi <- sei^2
  k <- length(yi)
  objective <- function(tau2) {
    wi <- 1 / (vi + tau2)
    mu <- sum(wi * yi) / sum(wi)
    sum(log(vi + tau2)) + log(sum(wi)) + sum(wi * (yi - mu)^2)
  }
  upper <- max(1, 10 * stats::var(yi))
  opt <- optimize(objective, interval = c(0, upper))
  tau2 <- opt$minimum
  wi <- 1 / (vi + tau2)
  mu <- sum(wi * yi) / sum(wi)
  se_mu <- sqrt(1 / sum(wi))
  q <- sum((1 / vi) * (yi - sum(yi / vi) / sum(1 / vi))^2)
  i2 <- if (q > 0) max(0, (q - (k - 1)) / q) * 100 else 0
  pred_se <- sqrt(tau2 + se_mu^2)
  tibble(
    countries = k, pooled_estimate = mu, pooled_SE = se_mu,
    CI_lo = mu - 1.96 * se_mu, CI_hi = mu + 1.96 * se_mu,
    tau2 = tau2, tau = sqrt(tau2), Q = q, Q_df = k - 1,
    Q_p = pchisq(q, df = k - 1, lower.tail = FALSE), I2_percent = i2,
    prediction_lo = mu - 1.96 * pred_se,
    prediction_hi = mu + 1.96 * pred_se
  )
}

meta_input <- country_results %>%
  filter(status %in% c("ok", "singular_random_slope"),
         is.finite(estimate), is.finite(SE), SE > 0)
meta_results <- meta_input %>%
  group_by(estimand) %>%
  group_modify(~reml_meta(.x$estimate, .x$SE)) %>%
  ungroup()
write.csv(meta_results, file.path(out_dir, "share_country_random_effects_summary.csv"), row.names = FALSE)

# Compare the prespecified CHARLS revised-M0 estimate with the distribution of
# SHARE country estimates. This is descriptive and not a causal comparison.
charls <- read.csv(file.path(out_dir, "revised_primary_combined_results.csv")) %>%
  filter(comparison == "China") %>%
  select(estimand, charls_estimate = estimate, charls_SE = SE)
distribution_compare <- meta_input %>%
  inner_join(charls, by = "estimand") %>%
  group_by(estimand) %>%
  summarise(
    charls_estimate = first(charls_estimate),
    countries = n(),
    country_min = min(estimate), country_median = median(estimate),
    country_max = max(estimate),
    countries_at_or_above_charls = sum(estimate >= first(charls_estimate)),
    .groups = "drop"
  )
write.csv(distribution_compare, file.path(out_dir, "charls_vs_share_country_distribution.csv"), row.names = FALSE)

plot_data <- meta_input %>%
  filter(estimand == "average_annual_change_difference_0_to_4") %>%
  arrange(estimate) %>%
  mutate(country = factor(country, levels = country))
charls_avg <- charls %>%
  filter(estimand == "average_annual_change_difference_0_to_4") %>%
  pull(charls_estimate)
p <- ggplot(plot_data, aes(estimate, country)) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_vline(xintercept = charls_avg, colour = "#B2182B", linetype = "dashed", linewidth = 0.7) +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), width = 0.18,
                colour = "grey45", orientation = "y") +
  geom_point(size = 2, colour = "#2166AC") +
  labs(
    x = "Mid versus low education difference in average annual memory change (points/year)",
    y = NULL,
    caption = "Country estimates use an identical model; dashed red line is the CHARLS revised-M0 estimate."
  ) +
  theme_classic(base_size = 10)
ggsave(file.path(out_dir, "share_country_forest_plot.png"), p, width = 7.2, height = 7.8, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "share_country_forest_plot.pdf"), p, width = 7.2, height = 7.8, device = cairo_pdf)

cat("\nCountry-model results\n"); print(as.data.frame(country_results), row.names = FALSE)
cat("\nRandom-effects synthesis\n"); print(meta_results)
cat("\nCHARLS versus SHARE-country distribution\n"); print(distribution_compare)
cat("\nSHARE country synthesis complete\n")
