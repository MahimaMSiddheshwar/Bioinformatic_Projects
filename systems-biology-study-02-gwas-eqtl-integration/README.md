# Systems Biology Study 02: GWAS and eQTL Integration for Alzheimer Disease

## Overview
This project integrates Alzheimer GWAS signals with eQTL evidence to prioritize significant variants, mapped genes, and genomic context.

## Objective
- Filter significant SNPs using genome-wide significance thresholds.
- Map significant loci to associated genes.
- Evaluate distance and positional context of SNP-gene relationships.
- Integrate BRAINEAC eQTL evidence and summarize associations.

## Data Inputs
- data/raw/gwas-association-file.tsv
- data/raw/EQTL.tsv

## Project Structure
- analysis/notebooks/gwas_eqtl_integration_analysis.ipynb - analysis workflow
- results/tables/MahimaMS_Alzheimers Significant SNPs.xlsx - prioritized SNPs
- results/tables/MahimaMS_eQTL Associated Gene.xlsx - eQTL-associated genes
- results/tables/eQTL Associated Gene Results.xlsx - additional export

## Requirements
- Python
- pandas, matplotlib

## How to Run
1. Install required Python libraries.
2. Open analysis/notebooks/gwas_eqtl_integration_analysis.ipynb.
3. Update data file paths if necessary.
4. Execute notebook cells in order.

## Outputs
- Filtered SNP table for Alzheimer-associated loci
- eQTL-linked gene summary tables
- Chromosome and distance-based visual summaries

## Notes
- This project focuses on integrative interpretation rather than model training.
