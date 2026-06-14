# Differential Expression Results Directory

This directory contains statistical output tables for the multivariable covariate-adjusted OLS regressions and downstream logistic classification modeling.

## Included Files
- `mirna_de_results.csv` / `mrna_de_results.csv`: OLS coefficients, p-values, and adjusted q-values for the simulated dataset.
- `real_mirna_de_results.csv` / `real_mrna_de_results.csv`: Multivariable adjusted OLS regression metrics for the GSE16759 cohort (controlling for Age, Sex, PMI, and RIN).
- `disease_classification_metrics.csv` / `disease_classification_metrics.json`: Performance characteristics (coefficients, p-values, ROC-AUC) of univariate penalized Logistic Regression classifiers evaluated on technical-adjusted residuals.

## Git Ignore Status
All result tables are generated dynamically at runtime and are ignored via `.gitignore`.
