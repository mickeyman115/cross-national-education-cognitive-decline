# Education and episodic-memory change across four ageing cohorts

This repository contains the locked analysis code and aggregate outputs for a
longitudinal comparison of education and episodic-memory change in
CHARLS, HRS, ELSA, and SHARE. The primary analysis included 144,642
participants and 542,426 episodic-memory observations.

## Scope

The primary estimand is the intermediate-minus-low education difference
in average annual episodic-memory change from cohort entry to year 4. The model
combines first-return and subsequent-response weights, demographic and retest
trajectory terms, and participant-specific random intercepts and time slopes.
The repository also implements a secondary year-4 instantaneous-rate estimand
and a common-specification synthesis across supported SHARE countries.

Scripts 04–10 preserve earlier socioeconomic and economic-position analyses for
audit history. They use an earlier exploratory specification not used for the
reported analyses and are not the basis for attenuation or mediation claims.

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
Rscript analysis/11_primary_without_first_return.R
Rscript analysis/12_first_return_selection.R
Rscript analysis/13_share_country_synthesis.R
Rscript analysis/14_final_outputs.R
Rscript analysis/15_make_manuscript_figures.R
```

Alternatively, `Rscript run_pipeline.R` executes the primary analysis
through the final aggregate result tables. Figure generation is kept separate.

The mixed models are computationally intensive. The aggregate locked outputs
are supplied in `results/` so readers can audit the reported estimates without
receiving restricted participant records.

## Software

The locked run used R 4.6.0. Package and platform details are recorded in
`results/final_sessionInfo.txt`. Principal packages include `haven`, `dplyr`,
`tidyr`, `mice`, `lme4`, `lmerTest`, `emmeans`, `ggplot2`, and `patchwork`.

## Interpretation boundary

This observational analysis does not establish a causal effect of education.
First-return and conditional-response weighting address selection associated
with measured variables under modelling assumptions; they do not prove absence
of selection or mortality bias. Study-specific education categories do not
establish equivalent schooling dose, quality, or social meaning across cohorts.

## Licence and citation

Code is released under the MIT License. Cohort data remain subject to their own
terms and are not covered by this licence. Please use `CITATION.cff` when citing
the software; the article citation will be added after publication.
