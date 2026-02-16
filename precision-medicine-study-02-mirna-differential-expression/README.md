# Precision Medicine Study 02: SARS-CoV-2 miRNA Differential Expression

## Overview
This project implements an RNA-seq workflow to analyze human respiratory cell responses to SARS-CoV-2 infection, from preprocessing to differential expression and GO enrichment.

## Objective
- Process raw RNA-seq data and quantify expression.
- Compare control vs infected and 24H vs 72H conditions.
- Identify differentially expressed gene and miRNA signals.
- Interpret results using GO enrichment and volcano plots.

## Study Design
- Conditions: Control vs SARS-CoV-2 infected
- Timepoints: 24H and 72H

## Project Structure
- analysis/scripts/mirna_differential_expression_analysis.R - DEG and enrichment workflow
- analysis/scripts/analysis_shell_commands.txt - command-line pipeline reference
- data/raw/miRNA_counts.txt - count matrix input
- results/tables/ - DEG and GO result tables
- results/figures/ - volcano and GO plots
- docs/report/mirna_differential_expression_report.docx - supporting write-up

## Pipeline Summary
1. Download SRA samples and convert to FASTQ.
2. Run quality control with FastQC and MultiQC.
3. Trim reads and align to hg38 with HISAT2.
4. Quantify with featureCounts to generate miRNA_counts.txt.
5. Perform DEG analysis with edgeR and limma.
6. Perform GO enrichment with clusterProfiler.

## Requirements
- Bash tools: sra-toolkit, fastqc, HISAT2, samtools, subread, multiqc
- R packages: edgeR, limma, clusterProfiler, org.Hs.eg.db, EnhancedVolcano

## Outputs
- results/tables/DEG_Comparison_Control_vs_Infected.csv
- results/tables/DEG_Comparison_24H_vs_72H.csv
- results/tables/GO_enrichment_Control_vs_Infected.csv
- results/tables/GO_enrichment_24H_vs_72H.csv
- results/figures/volcano_control_vs_infected.png
- results/figures/volcano_24h_vs_72h.png
- GO barplot outputs in results/figures/

## Notes
- Original command-level details are preserved in analysis/scripts/analysis_shell_commands.txt.
- Several commands contain environment-specific paths and should be adapted locally.
