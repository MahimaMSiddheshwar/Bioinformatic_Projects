# High-Throughput Study 02: RNA Decay and Half-Life Analysis

## Overview
This project estimates transcript half-life values from yeast time-course data across three replicates and performs functional enrichment for stable and unstable gene groups.

## Objective
- Fit decay models to estimate transcript half-lives.
- Compute average half-life across replicates.
- Identify top 10 percent and bottom 10 percent stability groups.
- Perform GO-based enrichment analysis.

## Data Inputs
- data/raw/decay_timecourse_raw_data.csv

## Method Summary
1. Clean missing values for each time course.
2. Fit exponential decay models per gene.
3. Compute replicate-specific half-lives and average half-life.
4. Rank genes by stability.
5. Run enrichment on top and bottom sets using gProfiler and GOTermFinder.

## Project Structure
- analysis/notebooks/rna_decay_half_life_analysis.ipynb - complete analysis workflow
- data/raw/ - original time-course dataset
- data/processed/ - cleaned and interpolated datasets
- results/tables/ - half-life and ranked gene outputs
- results/figures/ - enrichment result figures
- docs/report/rna_decay_stability_report.docx - supporting write-up

## Requirements
- Python
- pandas, numpy, scipy
- Jupyter Notebook

## How to Run
1. Install dependencies.
2. Open analysis/notebooks/rna_decay_half_life_analysis.ipynb.
3. Update data paths if needed.
4. Run all cells in sequence.

## Outputs
- results/tables/Average_Half_Lives_Results.csv
- results/tables/Top_10_Percent_Genes.csv
- results/tables/Bottom_10_Percent_Genes.csv
- replicate-level half-life tables in results/tables/
- enrichment figures in results/figures/

## Notes
- Some notebook paths are hardcoded and must be adjusted for local execution.
- Full code comments are included inside the notebook.
