# Figure Visualizations Directory

This directory holds the high-resolution vector diagrams (SVG format) illustrating differential expression, residual correlation stability, interactome bipartite networks, and protein-protein networks.

## Included Files
- `differential_expression_volcano.svg` / `real_volcano_plots.svg`: Volcano plots representing AD-associated transcript alterations.
- `confounder_adjusted_correlation.svg` / `real_adjusted_correlation.svg`: Residual correlation scatter plot for key target pairs.
- `real_interactome_network.svg`: Bipartite network connecting microRNAs to their target genes.
- `target_ppi_network.svg`: Protein-Protein interaction network visual retrieved programmatically from the STRING database.

## Git Ignore Status
All visual SVG assets are ignored via `.gitignore` to maintain a light-weight repository. They are plotted and saved automatically by the script engines.
