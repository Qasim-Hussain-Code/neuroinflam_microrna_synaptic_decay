# Walkthrough - Post-Doc miRNA-mRNA PFC Transcriptomics Pipeline

This guide outlines the components of the workspace and provides instructions for running the analysis pipeline locally.

---

## 1. Project Components Created

The following workspace files have been set up in `c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis`:

- **[environment.yml](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/environment.yml)**: Optimized Conda environment configuration pinning standard scientific packages (`numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, and `seaborn`) to ensure reproducible runs.
- **[target_database.json](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/target_database.json)**: A dictionary mapping brain-enriched microRNAs to their experimentally validated mRNA target transcripts (e.g., *hsa-miR-132* targeting *ITPKB*).
- **[run_mirna_pipeline.py](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/run_mirna_pipeline.py)**: The master Python engine running the analysis from end-to-end on simulated data.
- **[run_real_data_pipeline.py](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/run_real_data_pipeline.py)**: Real-world transcriptomics analysis script downloading, parsing, and running OLS and residual correlation models on GSE16759 dataset.
- **[.gitignore](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/.gitignore)**: Standard git exclusion criteria to keep environment folders and runtime-synthesized tables out of source control.
- **[README.md](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/README.md)**: Master documentation page detailing the research background, design matrix equations, and setup guides.
- **[LICENSE](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/LICENSE)**: Standard open-source MIT License terms.
- **[tests/test_pipeline.py](file:///c:/Users/User/Documents/bioinformatics_projects/miRNA_analysis/tests/test_pipeline.py)**: Unit test suite for validation of log-transform, design matrix dimensions, OLS parameters, and Fisher's exact test logic.

---

## 2. Analysis Phases Explained

When you run the script, it will execute ten distinct scientific phases:

### Phase 1: Cohort Synthesis & Biophysical Modeling (or Real Data Fetching)
- **Simulated Cohort ($N=120$)**: Generates a 120-donor cohort matching post-mortem characteristics of ROSMAP/Mayo studies, simulating differential autolysis degradation where mRNAs degrade rapidly ($\beta_{RIN} \approx 0.22$) while miRNAs associated with RISC complexes remain stable ($\beta_{RIN} \approx 0.05$).
- **Real-World Cohort ($N=8$ matched donors, GSE16759)**: Programmatically downloads parietal cortex microarray matrices, parses sample titles to pair matching donors, and extracts clinical metadata (`Age`, `Sex`, `PMI`, `Diagnosis`).

### Phase 2: Multivariable Covariate-Adjusted Regression (DE)
Fits a separate Ordinary Least Squares (OLS) model for each feature (miRNA and mRNA) controlling for technical autolysis and demographic confounders:
$$\log_2(\text{Expression} + 1) = \beta_0 + \beta_1(\text{AD\_Status}) + \beta_2(\text{RIN}) + \beta_3(\text{PMI}) + \beta_4(\text{Age}) + \beta_5(\text{Sex}) + \epsilon$$

*(For the real-world dataset, Age and PMI are centered around the cohort mean to prevent collinearity issues, and a Variance Inflation Factor (VIF) diagnostic is computed to verify covariate independence).*

### Phase 3: Confounder-Adjusted Residual Correlation with Bootstrap CIs
Instead of correlating raw expression, which yields false-positive signals due to co-degradation, it:
1. Regresses out all clinical and technical covariates.
2. Extracts the residuals ($\hat{\epsilon}$) representing variation independent of degradation or demographics.
3. Computes Pearson correlation on these residuals.
4. Runs **1,000-fold bootstrap resampling** to calculate 95% Confidence Intervals (CIs) to evaluate correlation stability under small cohort sizes.

### Phase 4: Over-Representation Analysis (ORA)
Performs a Fisher's Exact Test to determine if pairs showing significant negative correlation are statistically enriched for verified database targets from `target_database.json`.

### Phase 5: Pathway Analysis
Queries the Reactome Analysis Service API with the list of candidate targets to identify enriched physiological pathways (notably identifying the **Regulation of FOXO transcriptional activity** longevity axis).

### Phase 6: Figure Rendering
Generates high-resolution vector diagrams (SVG) in `results/figures/`:
1. `volcano_plots`: Volcano plots showing AD associations for miRNAs & mRNAs.
2. `adjusted_correlation`: Residual scatter plot showing negative correlation of a key validated target pair (e.g. *hsa-miR-9* vs. *SIRT1*).
3. `interactome_network`: Bipartite network diagram showing relationships between AD-dysregulated miRNAs (circular nodes) and target mRNAs (square nodes).

### Phase 7: STRING PPI Network Integration
Queries the STRING database API to find physical and functional interaction networks among target transcripts, verifying functional clustering and computing **PPI Enrichment p-value** (significant at $p = 0.0017$). Downloads the visual SVG network diagram from the STRING server.

### Phase 8: Diagnostic Disease Classification
Fits univariate **Logistic Regression models** on technical-adjusted residuals (which clear age, sex, and PMI technical confounders while retaining disease signals) to test if they can predict AD status, reporting **ROC-AUC scores** (achieving MAPK1 AUC = 0.938, hsa-miR-132 AUC = 0.875).

---

## 3. How to Run the Pipeline

Follow these commands in your terminal to initialize the environment and run the pipeline.

### Step 1: Initialize the Conda Environment
In your shell, run:
```bash
conda env create -f environment.yml
```

### Step 2: Activate the Environment
```bash
conda activate neuro_transcriptomics_env
```

### Step 3: Option A - Run the Real Data Analysis
```bash
python run_real_data_pipeline.py
```

### Step 4: Option B - Run the Simulated Analysis
```bash
python run_mirna_pipeline.py
```

---

## 4. Reviewing the Results

Once completed, you will find the following outputs in your workspace:

- **Raw Data**: Raw counts in `data/raw_expression/` and clinical variables in `data/clinical_metadata/`.
- **Differential Expression Metrics**: OLS regression parameters and adjusted Q-values under `results/differential_expression/`.
- **Classification Performance**: Univariate Logistic Regression ROC-AUC classification results in `results/differential_expression/disease_classification_metrics.csv` and `.json`.
- **Correlation Reports**: Pearson $r$ coefficients, bootstrap 95% CIs, robustness flags, and validation status under `results/correlation_analysis/`.
- **Target Enrichment statistics**: Contingency counts and Fisher's exact $p$-values in `results/correlation_analysis/target_enrichment_results.json`.
- **STRING PPI Networks**: STRING connection weights table in `results/correlation_analysis/target_ppi_network.tsv` and enrichment summary in `target_ppi_enrichment.json`.
- **Figures**: Clean vector graphics under `results/figures/`, including the STRING network diagram `target_ppi_network.svg`.
