# ======================================================================================
#                               Mahima M Siddheshwar - Assignment 03
# ======================================================================================


# 1. Install Required Packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute", force = TRUE)
BiocManager::install("preprocessCore", force = TRUE)
install.packages("WGCNA")




# ============ Note: Path to file are hardcoded please update the file path before running the code. =====================



# Disable automatic conversion of strings to factors (R < 4.0 compatibility)
options(stringsAsFactors = FALSE)

# 2. Load and Clean Data

# Load Required Packages
library(WGCNA)

expr_raw <- read.csv("C:/Users/mmsid/OneDrive/Documents/1. Mahima_IUI/1. Semester/1. 4th Sem/1. System Biologics/Assignment/Assignment 03/HW03_expression.csv")
expr_raw[,1] <- make.unique(as.character(expr_raw[,1]))
rownames(expr_raw) <- expr_raw[,1]
expr_raw <- expr_raw[,-1]
expr_raw <- as.data.frame(t(expr_raw))

traits_raw <- read.csv("C:/Users/mmsid/OneDrive/Documents/1. Mahima_IUI/1. Semester/1. 4th Sem/1. System Biologics/Assignment/Assignment 03/HW03_Traits.csv")
traits_raw <- traits_raw[!is.na(traits_raw[,1]) & traits_raw[,1] != "", ]
traits_raw[,1] <- make.unique(as.character(traits_raw[,1]))
rownames(traits_raw) <- traits_raw[,1]
traits_raw <- traits_raw[,-1]

good_samples <- complete.cases(expr_raw) & complete.cases(traits_raw)
expr_clean <- expr_raw[good_samples, ]
traits_clean <- traits_raw[good_samples, ]

write.csv(expr_clean, "CLEAN_expression.csv")
write.csv(traits_clean, "CLEAN_traits.csv")

# ----------------------------------------------------------------------------------------


# 3. Choose Soft Threshold (Question 1)

powers = c(1:20)
sft = pickSoftThreshold(expr_clean, powerVector = powers, verbose = 5)

png("Mahima_Q1_R2_Plots.png", width = 1000, height = 500)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Fit Index (R^2)",
     type = "n", main = "Scale-Free Topology Fit")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

softPower <- 13


# -------------------------------------------------------------------------------------


# 4. Gene Clustering via TOM (Question 2)

adjacency <- adjacency(expr_clean, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

png("Mahima_Q2_TOM_Dendrogram.png", width = 1000, height = 600)
plot(geneTree, xlab="", sub="", main = "Gene Clustering on TOM-based Dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# -----------------------------------------------------------------------------------------


# 5. Module Detection and Merging (Question 3)

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

png("Mahima_Q3_ModuleColors.png", width = 1000, height = 600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Merge modules
MEList <- moduleEigengenes(expr_clean, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

merge <- mergeCloseModules(expr_clean, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

png("Mahima_Q3_MergedModuleDendrogram.png", width = 1000, height = 600)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged Modules"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# ---------------------------------------------------------------------------------------


# 6. Calculate Module Eigengenes & Correlate with Traits (Question 4)

MEList <- moduleEigengenes(expr_clean, colors = mergedColors)
mergedMEs <- MEList$eigengenes  # This contains the eigengenes

# Correlate each module eigengene with traits
moduleTraitCor <- cor(mergedMEs, traits_clean, use = "p")

# Calculate p-values for significance of each correlation
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(expr_clean))

# Save correlation values and p-values
write.csv(moduleTraitCor, "Mahima_ModuleTrait_Correlations.csv")
write.csv(moduleTraitPvalue, "Mahima_ModuleTrait_PValues.csv")

# Calculate the dissimilarity of module eigengenes
MEDiss <- 1 - cor(mergedMEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Save the dendrogram 
png("Mahima_Q4_Module_Eigengene_Clustering.png", width = 1000, height = 600)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "", cex = 0.9)
abline(h = 0.25, col = "red")  
dev.off()


# ---------------------------------------------------------------------------------------


# 7. Module–Trait Association (Question 5)

moduleTraitCor <- cor(mergedMEs, traits_clean, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(expr_clean))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

png("Mahima_Q5_ModuleTraitHeatmap.png", width = 1200, height = 800, res = 150)

par(mar = c(6, 9, 4, 2))  

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits_clean),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,  
               cex.text = 0.6,
               zlim = c(-1,1),
               main = "Module–Trait Relationships")

dev.off()

#----------------------------------------------------------------------------------------------#



