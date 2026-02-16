# Systems Biology Study 04: Estrogen Response GSEA and Enrichment Mapping

## Overview
This project investigates estradiol treatment effects in MCF-7 breast cancer cells at 12h and 24h using differential expression, GSEA, Cytoscape enrichment maps, and heatmap visualization.

## Objective
- Identify pathways enriched after estrogen treatment.
- Compare pathway behavior across 12h and 24h timepoints.
- Visualize pathway relationships as enrichment maps.
- Summarize top core-enriched genes with heatmaps.

## Dataset
- GEO accession: GSE11352
- Platform: Affymetrix Human Genome U133A (GPL96)
- Samples: treated vs untreated at 12h and 24h

## Method Summary
1. Preprocess expression data and run limma differential analysis.
2. Build ranked gene lists for each timepoint.
3. Run GSEA using GO biological process gene sets.
4. Create enrichment maps in Cytoscape with EnrichmentMap and AutoAnnotate.
5. Generate heatmaps of top core-enriched genes.

## Project Structure
- analysis/scripts/gsea_enrichment_mapping_analysis.R - R workflow
- results/figures/ - dot plots, enrichment maps, heatmap
- docs/report/ - supporting documentation

## Tools and Packages
- R packages: limma, clusterProfiler, fgsea, org.Hs.eg.db, pheatmap, msigdbr
- Software: GSEA Desktop, Cytoscape 3.10.3, EnrichmentMap, AutoAnnotate

## Outputs
- results/figures/GSEA_12h_dotplot.png
- results/figures/GSEA_24h_dotplot.png
- results/figures/Enrichment_Map_12h.png
- results/figures/Enrichment_Map_24h.png
- results/figures/Top50_GSEA_Expression_Heatmap.png

## Key Interpretation
- Early response at 12h highlights immediate transcriptional and RNA-processing processes.
- Later response at 24h shows stronger enrichment in cell cycle and chromatin programs.
- Combined R and Cytoscape workflow provides both statistical and network-level interpretation.
