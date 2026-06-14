# Curated Expression Matrices Directory

This directory contains intermediate aligned, normalized, and log2-transformed expression matrices for both mRNA and miRNA transcripts.

## Included Files
- `mrna_raw_counts.tsv` / `mirna_raw_counts.tsv`: Aligned raw counts for target transcripts.
- `mrna_normalized_tpm.tsv` / `mirna_normalized_tpm.tsv`: Normalized (TPM-equivalent) expression matrices.
- `real_mrna_parietal_tpm.tsv` / `real_mirna_parietal_tpm.tsv`: Normalized parietal cortex expression matrices aligned for the matched 8 donors from GSE16759.

## Git Ignore Status
All curated matrices are ignored via `.gitignore` to maintain a lean repository. They are generated dynamically during the pipeline's preprocessing phase.
