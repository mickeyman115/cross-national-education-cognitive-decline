# Analytic-code release decision

## Locked decision

Public code deposition is desirable and should be completed before publication, but a fabricated or unsafe repository must not be supplied at initial submission. The individual-level CHARLS, HRS, ELSA, and SHARE data must never be placed in the repository.

## Initial-submission statement

Individual-level data are controlled by the respective cohort repositories and cannot be redistributed by the authors. Researchers can obtain CHARLS, HRS, ELSA, and SHARE data subject to their registration, application, and data-use conditions. Aggregate results supporting the findings are included in the Article and appendix. Analytic code will be prepared for public release without cohort data, and the repository URL will be added before publication.

## Release gate

Before making the repository public:

1. Retain only the canonical end-to-end analysis and figure scripts; exclude debug, pilot, and superseded scripts.
2. Remove absolute local paths, credentials, data-use tokens, temporary logs, and any row-level extracts.
3. Add a README with cohort acquisition routes, required data products, script order, expected outputs, R version, and package versions.
4. Add a synthetic schema or empty data dictionary sufficient to explain required variables without redistributing data.
5. Confirm that all exported examples, logs, and diagnostics contain no participant-level values.
6. Assign a stable public URL and release tag; preferably archive the tagged release with a DOI-capable repository.
7. Replace `[REPOSITORY URL TO BE ADDED]` in both manuscripts and submission forms.

## Claim boundary

Sharing code improves transparency but cannot by itself make the analysis fully executable for readers who have not independently obtained and harmonised the four restricted cohort datasets.
