# High-Throughput Study 03: Comparative Operon Prediction

## Overview
This project predicts operons across multiple organisms using adjacent co-directional gene rules and intergenic distance criteria.

## Objective
- Predict operon groups for E. coli, B. subtilis, Halobacterium, Synechocystis, and Hoatzin dataset.
- Export organism-specific operon prediction outputs.

## Data Inputs
- data/raw/E_coli_K12_MG1655.ptt
- data/raw/B_subtilis_168.ptt
- data/raw/Halobacterium_NRC1.ptt
- data/raw/Synechocystis_PCC6803_uid159873.ptt
- data/raw/Hoatzin.gff

## Method Summary
1. Load genome annotation files.
2. Sort genes by genomic position.
3. Group adjacent co-directional genes with short intergenic distances.
4. Export predicted operon sets per organism.

## Project Structure
- analysis/notebooks/comparative_operon_prediction.ipynb - main notebook
- data/raw/ - input annotation files
- results/tables/ - predicted operon outputs
- docs/report/operon_prediction_report.docx - supporting write-up

## Requirements
- Python
- pandas
- Jupyter Notebook

## How to Run
1. Install dependencies.
2. Open analysis/notebooks/comparative_operon_prediction.ipynb.
3. Update local data paths if required.
4. Execute notebook cells in order.

## Outputs
- results/tables/predicted_operons_B_subtilis.txt
- results/tables/predicted_operons_E_coli.txt
- results/tables/predicted_operons_Halobacterium.txt
- results/tables/predicted_operons_synechocystis.txt
- results/tables/predicted_operons_hoatzin.txt

## Notes
- Input paths may require updates before execution.
