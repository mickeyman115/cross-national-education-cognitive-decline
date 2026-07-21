###############################################################################
# Shared configuration for the public code release.
#
# Run scripts from the repository root. Cohort data are not distributed here.
# Set the four environment variables below to files obtained from the data
# custodians. EDUCOG_PROJECT_DIR defaults to the repository root.
###############################################################################

PROJECT_DIR <- normalizePath(
  Sys.getenv("EDUCOG_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_file <- function(variable) {
  value <- Sys.getenv(variable, unset = "")
  if (!nzchar(value) || !file.exists(value)) {
    stop(
      sprintf("Set %s to an authorised local cohort file before running this step.", variable),
      call. = FALSE
    )
  }
  normalizePath(value, mustWork = TRUE)
}

HRS_FILE <- function() required_file("HRS_FILE")
ELSA_FILE <- function() required_file("ELSA_FILE")
SHARE_FILE <- function() required_file("SHARE_FILE")
CHARLS_FILE <- function() required_file("CHARLS_FILE")

