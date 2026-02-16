# High-Throughput Study 01: Pearson Correlation Matrix Analysis

## Overview
This project computes Pearson correlation matrices from two expression datasets and visualizes correlation structure using heatmaps.

## Objective
- Compute Pearson correlation for Matrix 1 and Matrix 2.
- Compare correlation patterns between both matrices.
- Generate publication-ready heatmap outputs.

## Data Inputs
- data/raw/Mahima_Matrix1.txt
- data/raw/Mahima_Matrix2.txt

## Method Summary
1. Load tabular data into Python.
2. Compute Pearson correlation matrices using pandas and numpy.
3. Plot heatmaps with matplotlib and seaborn.
4. Export correlation tables and figures.

## Project Structure
- analysis/notebooks/matrix_pca_analysis.ipynb - main notebook workflow
- data/raw/ - source matrix files
- results/tables/ - computed correlation tables
- results/figures/ - heatmap visualizations
- docs/report/matrix_pca_report.docx - supporting write-up

## Requirements
- Python
- pandas, numpy, matplotlib, seaborn
- Jupyter Notebook

## How to Run
1. Install dependencies.
2. Open analysis/notebooks/matrix_pca_analysis.ipynb.
3. Update file paths if needed.
4. Run notebook cells in order.

## Outputs
- results/tables/PC_Matrix1_data.csv
- results/tables/PC_Matrix2_data.csv
- results/figures/PC_Matrix 1.png
- results/figures/PC_Matrix 2.png
- results/figures/PC_Matrix 1_Matrix 2_Combined.png

## Notes
- Some paths in notebook cells may be hardcoded; update to local paths before execution.
