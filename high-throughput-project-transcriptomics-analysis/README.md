# High-Throughput Project: Biomarker Discovery for Pregnancy-Associated Breast Cancer

## Overview
This project analyzes gene expression data to identify potential biomarkers linked to Pregnancy-Associated Breast Cancer (PABC) using transcriptomic differential expression and enrichment workflows.

## Objective
- Process and normalize expression data.
- Identify differentially expressed genes.
- Interpret gene lists with GO and KEGG pathway enrichment.
- Produce clear visual outputs for biological insight.

## Dataset
- GEO dataset: GSE31192

## Method Summary
1. Acquire expression data from GEO.
2. Perform preprocessing and normalization.
3. Run differential expression analysis.
4. Perform functional enrichment using GO and KEGG.
5. Generate figures and report outputs.

## Project Structure
- analysis/scripts/transcriptomics_analysis_pipeline.R - main R workflow
- results/figures/ - PCA, box plot, heatmaps, GO plots
- docs/report/Mahima MS_HTP_Project Report.pdf - project report

## Requirements
- R (tested with R 4.3.x)
- R packages: GEOquery, dplyr, ggplot2, pheatmap, ggrepel, limma, clusterProfiler, org.Hs.eg.db
- RStudio recommended

## How to Run
1. Install required R packages.
2. Open analysis/scripts/transcriptomics_analysis_pipeline.R in RStudio.
3. Update local paths and output directories.
4. Run script sections sequentially.

## Outputs
- Differential expression results
- Pathway enrichment summaries
- Figures in results/figures/ (PCA, heatmap, box plot, GO plots)

## Notes
- Internet access is required for GEO data retrieval via GEOquery.
- Update hardcoded file paths before running the script.
