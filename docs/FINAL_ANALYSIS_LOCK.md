# Final analysis lock

**Lock date:** 2026-07-21
**Status:** primary analysis complete
**Supersedes:** the earlier lock based on the 0.220 year-4 estimate

## Primary population and estimand

The baseline-eligible population contained 207,223 adults aged 50–100 years with complete age, sex, harmonised education, and at least one episodic-memory assessment. Longitudinal mixed-effects models used 144,642 repeated observers contributing 542,426 assessments. Stabilised first-return weights reweighted repeated observers towards the baseline-eligible population; these weights were multiplied by the existing conditional-response IPCW.

The primary estimand is the intermediate-minus-low education difference in average annual episodic-memory change from cohort entry to year 4. The year-4 instantaneous derivative is secondary.

## Primary results

| Comparison | Average annual difference, 0–4 years (95% CI), points/year |
|---|---:|
| CHARLS intermediate minus low | 0.093 (0.047 to 0.140) |
| HRS intermediate minus low | −0.019 (−0.047 to 0.008) |
| ELSA intermediate minus low | 0.005 (−0.022 to 0.032) |
| Pooled SHARE intermediate minus low | 0.004 (−0.007 to 0.014) |
| CHARLS minus HRS | 0.113 (0.059 to 0.167) |
| CHARLS minus ELSA | 0.088 (0.035 to 0.142) |
| CHARLS minus pooled SHARE | 0.090 (0.042 to 0.137) |

Holm-adjusted p values for the three cross-cohort contrasts were 0.00013, 0.00127, and 0.00049, respectively.

The secondary CHARLS year-4 instantaneous-rate contrast was 0.163 (95% CI 0.108 to 0.218). The corresponding CHARLS-minus-comparator contrasts were 0.165 for HRS, 0.148 for ELSA, and 0.180 for pooled SHARE.

## Selection analysis

Single-observation participants accounted for 23.4% of the CHARLS baseline-eligible population, 22.3% of HRS, 22.2% of ELSA, and 33.7% of SHARE. First-return models included education and sex in the numerator and flexible baseline age, flexible baseline memory, entry year, marital status, diabetes, hypertension, heart disease, and stroke in the denominator. The resulting first-return model c-statistics ranged from 0.660 to 0.893; cohort-specific effective sample-size ratios ranged from 0.925 to 0.968.

The model without first-return weighting estimated a CHARLS 0–4-year contrast of 0.09325; the combined-weight model estimated 0.09319. The combined-weight model converged, was non-singular, and had full fixed-effect rank (108/108). This supports robustness to selection explained by measured baseline variables. It does not eliminate unmeasured selection, pre-enrolment survival bias, or weight-estimation uncertainty.

## SHARE country heterogeneity

Country analyses were restricted to low and intermediate education and required at least 200 repeated observers, at least 50 people in each education category, at least 50 observations between years 3 and 5, and maximum follow-up of at least 4 years. Nineteen SHARE countries met these criteria. Every country used the same quadratic-time, flexible-age, sex-time, retest-education, random-intercept, and random-slope specification with the combined selection weight; no fallback model was used.

The random-effects pooled country estimate for average annual change was −0.021 points/year (95% CI −0.038 to −0.003; I²=38.9%; 95% prediction interval −0.072 to 0.030). Country estimates ranged from −0.134 to 0.042; none reached the CHARLS estimate of 0.093. Portugal had a boundary-singular random-slope fit and is retained with that diagnostic disclosed.

## Exploratory socioeconomic models

Earlier m=20 childhood-SEP and economic common-sample models remain archived as exploratory analyses. They were fitted under an earlier exploratory specification not used for the reported analyses and are not used to quantify attenuation of the primary estimand. They may be described only as hypothesis-generating supplementary evidence and must not be called mediation analyses.

## Claim ceiling

The locked conclusion is an observational cross-cohort association: the study-specific intermediate-versus-low education contrast in memory change was larger in CHARLS than in HRS, ELSA, pooled SHARE, and the supported SHARE-country distribution. The analysis does not demonstrate that education causes neuroprotection, prevents dementia or Alzheimer’s disease, establishes cognitive reserve, proves psychometric equivalence, or eliminates selection bias.

## Canonical outputs

- `results/final_revised_primary_results.csv`
- `output/first_return_selection_flow.csv`
- `output/first_return_selection_profile.csv`
- `output/first_return_weight_diagnostics.csv`
- `output/first_return_model_diagnostics.csv`
- `output/share_country_common_model_results.csv`
- `output/share_country_random_effects_summary.csv`
- `output/charls_vs_share_country_distribution.csv`

No further primary models, subgroups, or mechanism analyses are authorised for the initial submission unless a specific integrity check identifies a numerical or coding error.
