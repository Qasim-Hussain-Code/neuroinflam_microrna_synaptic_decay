# miRNA_analysis

[![Python Version](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![Conda](https://img.shields.io/badge/conda-enabled-green.svg)](https://conda.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An end-to-end post-doc tier transcriptomics framework investigating **miRNA-mRNA regulatory networks in Post-Mortem Alzheimer's Disease (AD) Prefrontal Cortex (PFC)**, designed to account for technical autolysis and tissue degradation.

---

## 1. Scientific Overview & Hypothesis

### Core Biological Question
> **"Does chronic microglial neuroinflammation in the human prefrontal cortex systematically correlate with a non-linear decay of post-synaptic density scaffold machinery ($DLG4$/$PSD\text{-}95$) when statistically adjusting for post-mortem degradation metrics?"**

### The Post-Mortem Confounder Shield
Human brain tissue collected post-mortem is highly volatile; technical covariates like post-mortem interval (PMI) and sample degradation (RNA Integrity Number, or RIN) skew genetic signals. To survive peer review, this framework applies a multi-variable OLS regression to isolate disease effects while neutralizing technical artifacts:

$$\log_2(\text{Expression}_{ij} + 1) = \beta_0 + \beta_1(\text{AD\_Status}_i) + \beta_2(\text{RIN}_i) + \beta_3(\text{PMI}_i) + \beta_4(\text{Age}_i) + \beta_5(\text{Sex}_i) + \epsilon_{ij}$$

Where:
* **`RIN` (RNA Integrity Number):** 1-10 scale measuring transcript degradation. (mRNAs show high sensitivity to degradation, whereas miRNAs associated with RISC complexes display higher technical resilience).
* **`PMI` (Post-Mortem Interval):** Time (hours) elapsed between death and flash-freezing.
* **`Age` & `Sex`:** Controls for demographic differences.

---

## 2. Directory Architecture

```text
miRNA_analysis/
├── environment.yml                  # Conda environment pins
├── .gitignore                       # Standard python and result ignore rules
├── target_database.json             # Curated brain miRNA targets database
├── run_mirna_pipeline.py            # Master simulation engine (with RIN/PMI covariates)
├── run_real_data_pipeline.py        # Real-world dataset engine (NCBI GEO GSE16759)
│
├── data/
│   ├── raw_expression/              # Raw count matrices (TSV format)
│   ├── clinical_metadata/           # Demographics and post-mortem covariates
│   └── curated/                     # Intermediate curated files
│
├── results/
│   ├── differential_expression/      # Covariate-adjusted OLS metrics and FDR Q-values
│   ├── correlation_analysis/         # Residual correlation coefficients & Fisher ORA
│   └── figures/                      # High-resolution vector diagrams (SVG format)
│
└── tests/
    └── test_pipeline.py             # Unit tests checking statistical engines
```

---

## 3. Quick Start Setup

### Step 1: Create the Environment
```bash
conda env create -f environment.yml
```

### Step 2: Activate the Environment
```bash
conda activate neuro_transcriptomics_env
```

### Step 3: Option A - Run the Real-World Dataset Analysis
This downloads the actual human brain microarray data from NCBI GEO GSE16759 (16 samples, 8 matched donors), maps probes, aligns samples, and runs OLS regressions and residual correlations on the real-world dataset:
```bash
python run_real_data_pipeline.py
```

### Step 3: Option B - Run the Simulated Covariate Analysis
This runs the high-sample-size (N=120) demographic simulation which explicitly models RIN and PMI decay kinetics:
```bash
python run_mirna_pipeline.py
```

### Step 4: Run the Tests
```bash
pytest tests/
```

---

## 4. Pipeline Processing Workflow

1. **Cohort Simulation & QC**: Generates a 120-donor cohort matching the technical and biological profile of ROSMAP/Mayo cohorts, modeling differential degradation rates for mRNA ($\beta \approx 0.22$) vs. miRNA ($\beta \approx 0.05$).
2. **OLS Differential Expression**: Regresses out clinical demographics (`Age`, `Sex`) and post-mortem confounders (`RIN`, `PMI`) to isolate the true AD effect ($\beta_1$). P-values are corrected using the Benjamini-Hochberg FDR procedure.
3. **Confounder-Adjusted Correlation**: Regresses out all design matrix covariates and calculates Pearson correlation on the OLS residuals ($\hat{\epsilon}$), verifying authentic biological co-expression.
4. **Target Enrichment Validation**: Matches negatively correlated miRNA-mRNA pairs against `target_database.json` and runs a Fisher's Exact Test to verify target enrichment.
5. **Vector Graphics Export**: Saves publication-grade volcano and residual correlation plots.
