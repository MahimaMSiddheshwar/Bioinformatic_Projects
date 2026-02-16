# Systems Biology Study 03: WGCNA Module-Trait Association

## Overview
This project applies Weighted Gene Co-Expression Network Analysis (WGCNA) to identify co-expression modules and evaluate module relationships with phenotypic traits.

## Objective
- Build co-expression modules from expression data.
- Select soft-threshold power using scale-free topology criteria.
- Correlate module eigengenes with traits.
- Summarize significant module-trait associations.

## Data Inputs
- data/processed/CLEAN_expression.csv
- data/processed/CLEAN_traits.csv

## Method Summary
1. Clean expression and trait matrices.
2. Select soft-thresholding power with diagnostic plots.
3. Build adjacency and TOM matrices.
4. Detect and merge modules.
5. Compute module eigengenes and module-trait correlations.
6. Visualize module-trait heatmaps and dendrograms.

## Project Structure
- analysis/scripts/wgcna_trait_association_analysis.R - main WGCNA script
- data/processed/ - cleaned expression and trait files
- results/tables/ - module-trait correlations and p-values
- results/figures/ - dendrogram and heatmap outputs

## Requirements
- R 4.2 or above
- WGCNA, BiocManager, impute, preprocessCore

## How to Run
1. Install required R packages.
2. Open analysis/scripts/wgcna_trait_association_analysis.R.
3. Set working directory and input paths.
4. Run script sections sequentially.

## Outputs
- results/tables/Mahima_ModuleTrait_Correlations.csv
- results/tables/Mahima_ModuleTrait_PValues.csv
- module dendrogram and module-trait heatmap figures

## Reference
- Method adaptation based on deneflab WGCNA tutorial sections 2.1 to 3.1.
