# Clinical Metadata Directory

This directory is designated for clinical and demographic characteristics of the donors in the analyzed cohorts.

## Included Files
- `clinical_metadata.tsv`: Simulated patient metadata (demographics, autolysis parameters) for the simulated cohort.
- `real_clinical_metadata.tsv`: Centered clinical metadata (Age, PMI, Sex, AD Status) for the matched GSE16759 cohort.

## Git Ignore Status
All tabular metadata files generated or parsed are ignored via `.gitignore` to prevent committing processed raw biological information. They are regenerated automatically during pipeline execution.
