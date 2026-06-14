# Correlation and Pathway Analysis Directory

This directory contains correlation coefficients, over-representation tests, pathway enrichments, and protein-protein interaction (PPI) tables.

## Included Files
- `mirna_mrna_correlation_results.csv` / `real_mirna_mrna_correlation_results.csv`: Pearson correlation coefficients, p-values, FDR q-values, and 95% bootstrap confidence intervals for expression residuals.
- `target_enrichment_results.json`: Hypergeometric ORA (Fisher's exact test) contingency tables and statistics evaluating microRNA target enrichment.
- `reactome_pathway_enrichment.csv`: Programmatic Reactome pathway enrichment outputs mapping target transcripts.
- `target_ppi_enrichment.json` / `target_ppi_network.tsv`: Functional and physical network connection details retrieved from the STRING database API.

## Git Ignore Status
All tabular and structured JSON analysis results are ignored via `.gitignore` and are populated dynamically when running the pipeline.
