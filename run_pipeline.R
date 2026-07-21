#!/usr/bin/env Rscript

# Complete locked workflow. Run from the repository root after setting the four
# authorised cohort-file environment variables documented in README.md.

steps <- c(
  "analysis/01_data_prep.R",
  "analysis/02_extract_covariates_ipcw.R",
  "analysis/03_build_ipcw.R",
  "analysis/04_extract_ses_wealth.R",
  "analysis/05_impute_common_sample.R",
  "analysis/06_fit_mechanism_models.R",
  "analysis/07_extract_fullsample_child_ses.R",
  "analysis/08_primary_analysis.R",
  "analysis/09_sensitivity_analysis.R",
  "analysis/10_export_mechanism_results.R"
)

for (step in steps) {
  cat("\n=== Running", step, "===\n")
  status <- system2("Rscript", step)
  if (!identical(status, 0L)) {
    stop(sprintf("Pipeline stopped because %s returned status %s", step, status))
  }
}

cat("\nLocked analysis pipeline completed.\n")

