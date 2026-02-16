# Mahima M Siddheshwar

# Assignment - 04


# ---------------------------------------------
# Install & Load Dependencies
# ---------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "msigdbr", "DOSE", "fgsea"))
install.packages(c("tidyverse", "pheatmap"))
BiocManager::install("enrichplot")



# Load libraries
library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(DOSE)
library(msigdbr)
library(tidyverse)
library(pheatmap)

# ---------------------------------------------
# Step 1: Load GSE11352 Expression Data
# ---------------------------------------------
gse <- getGEO("GSE11352", GSEMatrix = TRUE)[[1]]
expr <- exprs(gse)
pdata <- pData(gse)

# Keep only 12h and 24h samples
pdata <- pdata %>%
  filter(grepl("12hr|24hr", title)) %>%
  mutate(
    time = case_when(grepl("12hr", title) ~ "12h", grepl("24hr", title) ~ "24h"),
    treatment = case_when(grepl("E2 treated", title, ignore.case = TRUE) ~ "treated", grepl("Untreated control", title, ignore.case = TRUE) ~ "untreated"),
    group = paste(time, treatment, sep = "_")
  )

# Match Expression Matrix with Sample Info
rownames(pdata) <- pdata$geo_accession
expr_sub <- expr[, rownames(pdata)]
if (max(expr_sub) > 100) { expr_sub <- log2(expr_sub + 1) }

# Step 3: Design Matrix + DE Analysis
pdata$group <- factor(pdata$group)
design <- model.matrix(~ 0 + group, data = pdata)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(expr_sub, design)
contrast.matrix <- makeContrasts(
  E2_12h = group12h_treated - group12h_untreated,
  E2_24h = group24h_treated - group24h_untreated,
  levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Step 4: Map Probe IDs → Gene Symbols
res_12h <- topTable(fit2, coef = "E2_12h", number = Inf, sort.by = "t")
res_12h$PROBE <- rownames(res_12h)
gpl <- getGEO("GPL96", AnnotGPL = TRUE)
gpl_data <- Table(gpl)
annot <- gpl_data[, c("ID", "Gene symbol")]
colnames(annot) <- c("PROBE", "SYMBOL")
annot <- annot[annot$SYMBOL != "", ]
annot <- annot[!duplicated(annot$PROBE), ]
res_12h_annot <- merge(res_12h, annot, by = "PROBE")
gene_12h <- bitr(res_12h_annot$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
res_12h_annot <- merge(res_12h_annot, gene_12h, by = "SYMBOL")
rank_df <- res_12h_annot %>% group_by(ENTREZID) %>% slice_max(order_by = abs(t), n = 1) %>% ungroup()
gene_list_12h <- rank_df$t
names(gene_list_12h) <- rank_df$ENTREZID
gene_list_12h <- sort(gene_list_12h, decreasing = TRUE)

# Load Gene Sets and Run GSEA for 12h
msigdb <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gmt <- msigdb %>% dplyr::select(gs_name, entrez_gene)
gsea_12h <- GSEA(geneList = gene_list_12h, TERM2GENE = gmt, pvalueCutoff = 0.05)
png("GSEA_12h_dotplot.png", width = 10, height = 8, units = "in", res = 300)
dotplot(gsea_12h, showCategory = 20, title = "GSEA Enrichment at 12h") +
  scale_color_continuous(low = "lightpink", high = "purple4") +
  theme_minimal()
dev.off()

# Step 5: GSEA for 24h
res_24h <- topTable(fit2, coef = "E2_24h", number = Inf, sort.by = "t")
res_24h$PROBE <- rownames(res_24h)
res_24h_annot <- merge(res_24h, annot, by = "PROBE")
gene_24h <- bitr(res_24h_annot$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
res_24h_annot <- merge(res_24h_annot, gene_24h, by = "SYMBOL")
rank_df_24h <- res_24h_annot %>% group_by(ENTREZID) %>% slice_max(order_by = abs(t), n = 1) %>% ungroup()
gene_list_24h <- rank_df_24h$t
names(gene_list_24h) <- rank_df_24h$ENTREZID
gene_list_24h <- sort(gene_list_24h, decreasing = TRUE)
gsea_24h <- GSEA(geneList = gene_list_24h, TERM2GENE = gmt, pvalueCutoff = 0.05)
png("GSEA_24h_dotplot.png", width = 10, height = 8, units = "in", res = 300)
dotplot(gsea_24h, showCategory = 20, title = "GSEA Enrichment at 24h") +
  scale_color_continuous(low = "lightpink", high = "purple4") +
  theme_minimal()
dev.off()

write.table(gsea_12h@result, file = "GSEA_12h_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gsea_24h@result, file = "GSEA_24h_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



# ------------------Enrichment Plot-----------------------

# Enrichplot library
library(enrichplot)
library(ggplot2)

# Subset top pathways to reduce clutter
gsea_12h_top <- gsea_12h
gsea_12h_top@result <- gsea_12h@result %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::slice_head(n = 30)


# Ensure pairwise similarity is computed
gsea_12h <- pairwise_termsim(gsea_12h)

# Save PNG to device
png("Enrichment_Map_12h.png", width = 12, height = 10, units = "in", res = 300)

emapplot(
  gsea_12h,
  showCategory = 20,       
  color = "p.adjust",       
  layout = "kk"             
)

dev.off()


# Compute similarity between terms
gsea_24h <- pairwise_termsim(gsea_24h)

# Save high-resolution PNG
png("Enrichment_Map_24h.png", width = 12, height = 10, units = "in", res = 300)

emapplot(
  gsea_24h,
  showCategory = 20,        
  color = "p.adjust",        
  layout = "kk"              
)

dev.off()



# ------------ Heatmap for Top 50 Core Enriched Genes ---------------

top_genes_12h <- gsea_12h@result %>% dplyr::filter(p.adjust < 0.25) %>% head(50) %>% pull(core_enrichment) %>% strsplit("/") %>% unlist()
top_genes_24h <- gsea_24h@result %>% dplyr::filter(p.adjust < 0.25) %>% head(50) %>% pull(core_enrichment) %>% strsplit("/") %>% unlist()
heatmap_genes <- unique(c(top_genes_12h, top_genes_24h))
heatmap_genes <- head(heatmap_genes, 50)
heatmap_gene_symbols <- bitr(heatmap_genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
heatmap_gene_symbols <- unique(na.omit(heatmap_gene_symbols))

# Build expr_annot from probe data
expr_annot <- expr_sub[rownames(expr_sub) %in% annot$PROBE, ]
rownames(expr_annot) <- annot$SYMBOL[match(rownames(expr_annot), annot$PROBE)]
expr_annot <- expr_annot[!duplicated(rownames(expr_annot)), ]

# Create annotation_col for pheatmap
annotation_col <- pdata %>% dplyr::select(time, treatment)
rownames(annotation_col) <- rownames(pdata)

# Subset expression matrix
expr_heatmap <- expr_annot[rownames(expr_annot) %in% heatmap_gene_symbols, ]
expr_heatmap <- expr_heatmap[rowSums(is.na(expr_heatmap)) == 0, ]
expr_heatmap <- expr_heatmap[rowSums(expr_heatmap) != 0, ]

# Plot heatmap
png("Top50_GSEA_Expression_Heatmap.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(expr_heatmap,
         scale = "row",
         annotation_col = annotation_col,
         fontsize_row = 6,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Heatmap for Treated vs Untreated at 12h & 24h",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()


#-----------------------------------------------------------------------------------------#