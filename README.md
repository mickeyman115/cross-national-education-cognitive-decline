# Education and cognitive decline across four ageing cohorts

This repository contains the locked analysis code and aggregate outputs for a
longitudinal comparison of the education-cognitive decline association in
CHARLS, HRS, ELSA, and SHARE. The primary analysis included 144,642
participants and 542,426 cognitive observations.

## Scope

The code estimates cohort-specific differences in the instantaneous cognitive
trajectory at four years between intermediate and low education groups using
quadratic-time linear mixed models. It also implements childhood socioeconomic
adjustment with multiple imputation, economic-variable mechanism exploration,
IPCW sensitivity analyses, CHARLS entry-wave analyses, and the manuscript
figures.

This is a code-transparent but data-restricted release. Individual-level data,
imputation objects, and fitted model objects are intentionally excluded.

## Repository structure

- `analysis/`: ordered R scripts and shared configuration
- `results/`: aggregate locked results and diagnostics
- `data/`: data-access instructions and a variable dictionary
- `docs/`: final analysis lock and public-release boundary

## Data preparation

Obtain the four harmonised cohort products listed in `data/README.md`. From the
repository root, set their authorised local paths:

```sh
export HRS_FILE="/authorised/path/randhrs1992_2020v2.dta"
export ELSA_FILE="/authorised/path/h_elsa_g3.dta"
export SHARE_FILE="/authorised/path/H_SHARE_f2.dta"
export CHARLS_FILE="/authorised/path/H_CHARLS_D_Data.dta"
```

Do not commit these files. The `.gitignore` excludes common restricted-data and
derived person-level formats.

## Analysis order

Run from the repository root:

```sh
Rscript analysis/01_data_prep.R
Rscript analysis/02_extract_covariates_ipcw.R
Rscript analysis/03_build_ipcw.R
Rscript analysis/04_extract_ses_wealth.R
Rscript analysis/05_impute_common_sample.R
Rscript analysis/06_fit_mechanism_models.R
Rscript analysis/07_extract_fullsample_child_ses.R
Rscript analysis/08_primary_analysis.R
Rscript analysis/09_sensitivity_analysis.R
Rscript analysis/10_export_mechanism_results.R
Rscript analysis/11_make_submission_figures.R
```

Alternatively, `Rscript run_pipeline.R` executes analysis steps 1–10 in order.
Figure generation is kept separate because it requires the fitted primary
model produced by the pipeline.

The complete pipeline is computationally intensive. The multiple-imputation
mixed models can require many hours on a standard workstation. The aggregate
locked outputs are supplied in `results/` so readers can audit the reported
estimates without receiving restricted participant records.

## Software

The locked run used R 4.6.0. Package and platform details are recorded in
`results/final_sessionInfo.txt`. Principal packages include `haven`, `dplyr`,
`tidyr`, `mice`, `lme4`, `lmerTest`, `emmeans`, `ggplot2`, and `patchwork`.

## Interpretation boundary

This observational analysis does not establish a causal effect of education.
The economic-variable models are mechanism-oriented sequential-adjustment
analyses, not formal mediation models. IPCW addresses observed predictors of
subsequent response among survivors and does not prove the absence of all
attrition or mortality bias.

## Licence and citation

Code is released under the MIT License. Cohort data remain subject to their own
terms and are not covered by this licence. Please use `CITATION.cff` when citing
the software; the article citation will be added after publication.
