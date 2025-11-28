# ==============================================================================
# EGPA scRNA-seq Analysis: Figure 2
# ==============================================================================
# Description: This script generates Figure 2A, 2B, 2C, 2F
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic
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

# Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
granu <- readRDS("./data/Fig2/granu.RDS")

# ==============================================================================
# FIGURE 2A: UMAP Plot of Granulocyte Clusters
# ==============================================================================
mypalette <- c("#A6CEE3", "#FB9A99", "#1F78B4", "#E31A1C", "#CAB2D6", "darkred",
               "#33A02C", "#FDBF6F", "darkgrey", "#7570B3", "#D95F02", "#6A3D9A",
               "#B15928", "#1B9E77", "#FF7F00", "#FFFF99", "black", "navy", "firebrick")
p <- DimPlot(granu, reduction = "umap", cols = mypalette, pt.size = 0.1,
             label = FALSE, raster = FALSE)

# ==============================================================================
# FIGURE 2B: Cell Proportions in Clusters
# ==============================================================================
pdata <- data.frame(table(granu$subject_id, granu$source, granu@active.ident))
pdata$TotalNumOfCelltype <- rep(unname(table(granu@active.ident)), each = 112)
a <- unique(granu@meta.data[, c("subject_id", "group")])
a <- a[order(a$subject_id), ]
pdata$group <- rep(a[, "group"], 56)
pdata$Proportion <- pdata$Freq / pdata$TotalNumOfCelltype
pdata$x <- paste(pdata$Var2, pdata$group, sep = "-")
pdata$x <- factor(pdata$x, levels = c("PBMC-C", "PBMC-A", "PBMC-E",
                                       "Biopsy-C", "Biopsy-A", "Biopsy-E",
                                       "BALF-C", "BALF-A", "BALF-E",
                                       "BALF_EOS-C", "BALF_EOS-A"))

ggplot(pdata, aes(x = Var3, y = Proportion, fill = x)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
  )

# ==============================================================================
# FIGURE 2C: Feature Plots of Gene Expression
# ==============================================================================
DefaultAssay(granu) <- "SCT"
genes <- c("TNF", "IL1B", "NLRP3", "CXCL2", "CXCL3",
           "IRF7", "IFIT1", "IFIT2", "ISG15", "MX1")
base_colors <- rev(brewer.pal(11, "RdYlBu")[1:7])
gradient_100 <- colorRampPalette(base_colors)(100)

FeaturePlot(granu, features = c("CXCL3"), slot = "data", cols = gradient_100,
            order = FALSE, min.cutoff = 0, max.cutoff = 0.75, pt.size = 0.1,
            reduction = "umap", raster = FALSE)

FeaturePlot(granu, features = genes, cols = rev(brewer.pal(11, "RdYlBu")[1:7]),
            reduction = "umap", min.cutoff = 0, max.cutoff = 5, pt.size = 1)

FeaturePlot(granu, features = genes, slot = "data", cols = base_colors,
            order = TRUE, min.cutoff = 0, max.cutoff = 5, pt.size = 2,
            reduction = "umap", ncol = 5)

# ==============================================================================
# FIGURE 2F: Dot Plot of cirEos and Eos DEGs
# ==============================================================================
granu@graphs <- list()
spe <- subset(granu, subset_2 == "cirEos")

head(granu@meta.data)
unique(granu@meta.data$group)

Idents(spe) <- "group"
DEGs <- FindMarkers(spe, ident.1 = "E", ident.2 = "C", logfc.threshold = 0)
RPgenes <- grep(pattern = "^RP", x = rownames(DEGs), value = TRUE)
MTgenes <- grep(pattern = "^MT-", x = rownames(DEGs), value = TRUE)
DEGs <- DEGs[which(!rownames(DEGs) %in% c(RPgenes, MTgenes)), ]

spe <- subset(granu, idents = c("cirEos", "Eos"))
spe$group_cluster <- paste(spe$group, spe$subset_2, sep = "_")
Idents(spe) <- "group_cluster"
spe <- subset(spe, idents = c("C_cirEos", "A_cirEos", "E_cirEos",
                               "C_Eos", "A_Eos", "E_Eos"))
spe@active.ident <- factor(spe@active.ident,
                           levels = rev(c("C_cirEos", "A_cirEos", "E_cirEos",
                                         "C_Eos", "A_Eos", "E_Eos")))

features <- c("CLC", "IL5RA", "IL17RA", "SELL", "RNASE2", "HIF1A", "SMPD3",
              "CASP8", "MAPK3", "PRSS33", "CD69", "ENO1", "APOE", "APOC1",
              "PLAUR", "SIGLEC8", "CCR3", "IL1RL1", "IL3RA", "IL4R", "IL13RA1")

DotPlot(spe, assay = "RNA", features = features, col.min = -1, col.max = 1,
        dot.min = 0, dot.scale = 8, scale.min = 0, scale.max = 50) +
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
  )

# ==============================================================================
# End of Script
# ==============================================================================
