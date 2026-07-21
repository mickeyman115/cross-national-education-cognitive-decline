#!/usr/bin/env Rscript

# Revised primary workflow. Run from the repository root after setting the four
# authorised cohort-file environment variables documented in README.md.

steps <- c(
  "analysis/01_data_prep.R",
  "analysis/02_extract_covariates_ipcw.R",
  "analysis/03_build_ipcw.R",
  "analysis/11_revised_primary_without_first_return.R",
  "analysis/12_first_return_selection.R",
  "analysis/13_share_country_synthesis.R",
  "analysis/14_final_revised_outputs.R"
)

for (step in steps) {
  cat("\n=== Running", step, "===\n")
  status <- system2("Rscript", step)
  if (!identical(status, 0L)) {
    stop(sprintf("Pipeline stopped because %s returned status %s", step, status))
  }
}

cat("\nRevised primary analysis pipeline completed.\n")
