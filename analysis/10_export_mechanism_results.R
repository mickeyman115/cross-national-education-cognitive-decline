library(dplyr)
source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
output_dir <- OUTPUT_DIR
proj_dir <- PROJECT_DIR

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

  mechanism_results <- pooled %>%
    mutate(
      is_cross = grepl("^China - ", study),
      p_value = ifelse(!is.na(p_value), p_value, 2 * pnorm(abs(Q_bar / SE_pooled), lower.tail = FALSE))
    ) %>%
    group_by(model) %>%
    mutate(p_holm = ifelse(is_cross, p.adjust(p_value[is_cross], method = "holm")[cumsum(is_cross)], NA_real_)) %>%
    ungroup()

  write.csv(mechanism_results, file.path(output_dir, "final_mechanism_results.csv"), row.names = FALSE)

  sens_summary <- pooled %>%
    filter(study == "China") %>%
    select(model, Q_bar, SE_pooled) %>%
    mutate(across(c(Q_bar, SE_pooled), ~ round(., 4)))
  write.csv(sens_summary, file.path(output_dir, "final_sensitivity_results.csv"), row.names = FALSE)
}
