# Raw Expression Data Directory

This directory is designated for raw expression matrices downloaded programmatically from the NCBI GEO database. 

## Included Files
- `GPL570_brief.txt`: Probe-to-gene mapping metadata for the GPL570 (mRNA) microarray platform.
- `GPL8757_brief.txt`: Probe-to-gene mapping metadata for the GPL8757 (miRNA) microarray platform.
- `mrna_series_matrix.txt.gz`: Gzipped series matrix file containing raw expression data for mRNA samples (GSE16759).
- `mirna_series_matrix.txt.gz`: Gzipped series matrix file containing raw expression data for miRNA samples (GSE16759).

## Git Ignore Status
To avoid committing large binary and text files to GitHub, all data matrix files in this directory are ignored via `.gitignore`. Running the pipeline scripts will automatically download and cache these files locally.
