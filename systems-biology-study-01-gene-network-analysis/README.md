# Systems Biology Study 01: Human Protein Interaction Network Analysis

## Overview
This project analyzes the human protein-protein interaction network using BioGRID data to evaluate topology, hub proteins, shortest paths, and centrality patterns.

## Objective
- Build the human PPI graph.
- Identify hub proteins from node degree.
- Evaluate scale-free behavior of network degree distribution.
- Assess shortest-path structure and centrality metrics.

## Dataset
- Source: BioGRID
- Input file: BIOGRID-ORGANISM-Homo_sapiens-4.4.218.tab3.txt

## Method Summary
1. Clean and preprocess interaction pairs.
2. Build an undirected graph in igraph.
3. Compute node degree and top hub proteins.
4. Fit degree distribution using a power-law model.
5. Compute shortest paths and centrality measures.

## Project Structure
- analysis/notebooks/gene_network_analysis.ipynb - main analysis notebook
- docs/report/gene_network_analysis_report.pdf - supporting report

## Requirements
- Python
- pandas, igraph, matplotlib, seaborn, powerlaw, numpy

## How to Run
1. Install dependencies.
2. Open analysis/notebooks/gene_network_analysis.ipynb.
3. Update input path to local BioGRID file.
4. Run cells sequentially.

## Outputs
- Network summary statistics
- Top hub node table
- Degree distribution and power-law plots
- Shortest path distribution plot
- Centrality summaries

## Results Summary
- The network shows scale-free behavior (gamma approximately 3.367).
- The top hub node has very high connectivity.
- Most shortest paths are concentrated in low path lengths.
