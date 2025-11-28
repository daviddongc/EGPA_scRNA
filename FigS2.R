# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S2
# ==============================================================================
# Description: This script generates Supplementary Figure 2E-G
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic 1
# granulomatosis with polyangiitis> publication
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(sctransform)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(purrr)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(rstatix)
library(xlsx)

# Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# Funtions
# ------------------------------------------------------------------------------
plotting_heatmap<-function(expr_mat){

m=expr_mat
m=m[!apply(m,1,sum)==0,]

#Row-center the data.
m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m=m[is.na(row.names(m)) == FALSE,]
m[is.nan(m)] = 0
library(Seurat)
m<-MinMax(m,-2,2)

#obtaining the order of peak expression of each gene:
name<-apply(m,1,function(x){
names(which.max(x))
})
order<-intersect(colnames(m),name)
row_names<-c()
for (i in 1:length(order)){
row_names<-c(row_names,names(which(name==order[i])))
}
m<-m[row_names,]

heatmap_matrix<-m

row_dist<-as.dist((1-cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)]<-1
bks<-seq(-2.1,2.1, by = 0.1)
library(RColorBrewer)
hmcols<-colorRampPalette(rev(c(brewer.pal(7,"RdBu"))))(length(bks)-1)

ph<-pheatmap(heatmap_matrix,useRaster=T,cluster_cols=FALSE,cluster_rows=FALSE,clustering_distance_rows=row_dist,
breaks=bks,color=hmcols)
return(ph)
}

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
combined <- readRDS("./data/FigS2/combined.RDS")

# ------------------------------------------------------------------------------
# Cell Type Annotation (Coarse Clustering)
# ------------------------------------------------------------------------------
mtd <- combined@meta.data
mtd <- within(mtd, {
  coarse_cluster.r1 <- as.character(seurat_clusters)
  coarse_cluster.r1[seurat_clusters %in% c(12, 39, 56)] <- "B cells"
  coarse_cluster.r1[seurat_clusters %in% c(50, 4, 31, 33, 5, 49, 22, 6, 11, 51)] <- "T & NK cells"
  coarse_cluster.r1[seurat_clusters %in% c(15, 30, 36, 54, 25, 57, 29, 8, 0, 45, 47)] <- "Macrophage"
  coarse_cluster.r1[seurat_clusters %in% c(27)] <- "Fibroblast"
  coarse_cluster.r1[seurat_clusters %in% c(35)] <- "Smooth muscle"
  coarse_cluster.r1[seurat_clusters %in% c(14, 19, 24, 42, 52, 46, 55, 10, 7, 16)] <- "Epithelial"
  coarse_cluster.r1[seurat_clusters %in% c(26, 53)] <- "Endothelial"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(41)] <- "Megakaryocyte"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(3, 17, 38)] <- "Neutrophils-1"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(43)] <- "Myeloblast"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(20)] <- "Mono DCs"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(18)] <- "Mast cells"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(40)] <- "pDC"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(34)] <- "Neuroendocrine"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(2, 9, 21, 23, 28, 37, 44)] <- "Eosinophils"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(1, 13, 48)] <- "Neutrophils-2"
  coarse_cluster.r1[coarse_cluster.r1 %in% c(32)] <- "Erythroid cells"
})

combined@meta.data <- mtd
Idents(combined) <- "coarse_cluster.r1"

# Set cell type order
combined@active.ident <- factor(combined@active.ident,
                                levels = c("Epithelial", "Neuroendocrine", "Fibroblast",
                                          "Smooth muscle", "Endothelial", "T & NK cells",
                                          "B cells", "Myeloblast", "Macrophage", "Mono DCs",
                                          "pDC", "Neutrophils-1", "Neutrophils-2",
                                          "Eosinophils", "Mast cells", "Megakaryocyte",
                                          "Erythroid cells"))

# Define color palette for cell types
mypalette <- c("#A6CEE3", "#FB9A99", "#1F78B4", "#E31A1C", "#CAB2D6", "darkred",
               "#33A02C", "#FDBF6F", "darkgrey", "#7570B3", "#D95F02", "#6A3D9A",
               "#B15928", "#1B9E77", "#FF7F00", "navy", "black", "firebrick", "#FFFF99")

# ==============================================================================
# FIGURE S2E: UMAP Plots by Source, Group, and Subject
# ==============================================================================

# UMAP by source
p <- DimPlot(combined, reduction = "umap", cols = c(brewer.pal(4, "Set2")),
            pt.size = 0.1, group.by = "source")
p
# UMAP by group
combined@meta.data$group <- factor(combined@meta.data$group, levels = c("C", "A", "E"))
p <- DimPlot(combined, reduction = "umap", cols = c("#5BBB93", "#377EB8", "#BA2833"),
            pt.size = 0.1, group.by = "group")
p
# UMAP by subject
p <- DimPlot(combined, reduction = "umap",
            cols = c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")),
            pt.size = 0.1, group.by = "subject_id")
p
# ==============================================================================
# FIGURE S2F: Gene Expression Correlation Analysis
# ==============================================================================

# Create subject-source identifier
combined@meta.data$subject_source <- paste(combined$subject_id, combined$source, sep = "_")
Idents(combined) <- "subject_source"

# Calculate average expression
avg <- AverageExpression(combined)

# Reset identity to cell types
Idents(combined) <- "coarse_cluster.r1"

# Prepare proportion data
pdata <- data.frame(table(combined$subject_id, combined$source, combined@active.ident))
pdata$Group <- rep(0, length(pdata[, 1]))
for (i in 1:length(pdata[, 1])) {
  pdata[i, "Group"] <- as.character(unique(combined@meta.data[which(combined$subject_id == pdata[i, 1]), "group"]))
}
pdata$col <- paste(pdata$Var1, pdata$Var2, sep = "_")

# Calculate correlation for each tissue and group combination
corExp <- data.frame()
for (i in c("PBMC", "Biopsy", "BALF", "BALF_EOS")) {
  for (j in c("C", "E", "A")) {
    col <- unique(pdata[which(pdata$Var2 == i & pdata$Group == j), "col"])
    col <- intersect(col, colnames(avg$SCT))
    mat <- avg$SCT[, col]
    cor <- unique(c(cor(mat)))
    cor <- cor[which(cor != 1)]
    cor <- data.frame(cbind(rep(i, length(cor)), rep(j, length(cor)), cor))
    corExp <- rbind(cor, corExp)
  }
}

# Plot correlation results
corExp$V1 <- factor(corExp$V1, levels = c("PBMC", "Biopsy", "BALF", "BALF_EOS"))
corExp$V2 <- factor(corExp$V2, levels = c("C", "A", "E"))
ggplot(corExp, aes(x = V1, y = as.numeric(as.character(cor)), fill = V2)) +
  geom_boxplot(outlier.color = NA) +
  ylim(0, 1) +
  theme_bw()

# ==============================================================================
# FIGURE S2G: Heatmap of Differentially Expressed Genes
# ==============================================================================

# Find all markers
FindAllMarkers(combined) -> combined.allmarkers

# Filter significant markers
a <- combined.allmarkers[which(combined.allmarkers[, 2] >= 0.5 & combined.allmarkers[, 5] <= 0.05), ]

# Select top 20 markers per cluster
top20 <- a %>%
  group_by(cluster) %>%
  dplyr::slice(1:20)

# Create heatmap
p <- DoHeatmap(combined, features = top20$gene, disp.min = -2, disp.max = 2, label = FALSE)
bks <- seq(-2.1, 2.1, by = 0.1)
heatmapColors <- colorRampPalette(rev(c(brewer.pal(7, "RdBu"))))(length(bks) - 1)
p + scale_fill_gradientn(colors = heatmapColors)

# Save results
write.csv(a, file = "./data/FigS2/differentially_expressed_genes.csv", row.names = FALSE)

# ==============================================================================
# End of Script
# ==============================================================================
