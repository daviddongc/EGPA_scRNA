# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S7
# ==============================================================================
# Description: This script generates Supplementary Figure S7
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

# ==============================================================================
# FIGURE S7: T/NK Cell and B Cell Receptor-Ligand Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# T/NK Cell Analysis
# ------------------------------------------------------------------------------
# Load T/NK cell data
TNK <- readRDS("./data/FigS7/S7_TNK.RDS")

# Define color palette
mypalette <- c("#A6CEE3", "#FB9A99", "#1F78B4", "#E31A1C", "#CAB2D6", "darkred",
               "#33A02C", "#FDBF6F", "darkgrey", "#7570B3", "#D95F02", "#6A3D9A",
               "#B15928", "#1B9E77", "#FF7F00", "#FFFF99", "black", "navy",
               "firebrick", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
               "#FFD92F", "#E5C494")

# S7_a: UMAP plot
p <- DimPlot(TNK, reduction = "umap", cols = mypalette, pt.size = 2, label = F)
ggsave("./data/FigS7/S7_a.pdf", p, width = 10, height = 10)

# S7_b: Cell proportions
# Calculate proportions of different individuals in each subtype
pdata <- data.frame(table(TNK$subject_id, TNK$source, TNK@active.ident))
pdata$TotalNumOfCelltype <- rep(unname(table(TNK@active.ident)), each = 112)
a <- unique(TNK@meta.data[, c("subject_id", "group")])
a <- a[order(a$subject_id), ]
pdata$group <- rep(a[, "group"], 88)
pdata$Proportion <- pdata$Freq / pdata$TotalNumOfCelltype
pdata$x <- paste(pdata$Var2, pdata$group, sep = "-")
pdata$x <- factor(pdata$x, levels = c("PBMC-C", "PBMC-A", "PBMC-E",
                                       "Biopsy-C", "Biopsy-A", "Biopsy-E",
                                       "BALF-C", "BALF-A", "BALF-E",
                                       "BALF_EOS-C", "BALF_EOS-A"))
p <- ggplot(pdata, aes(x = Var3, y = Proportion, fill = x)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()
ggsave("./data/FigS7/S7_b.pdf", p, width = 10, height = 10, dpi = 300)

# T cells Ligands and Receptors 1: Co-stimulatory and co-inhibitory molecules
genes <- c("CD40LG", "CD28", "CTLA4", "CD27", "HAVCR2", "ICOS", "PDCD1", "BTLA",
           "TNFRSF9", "CD5", "TNFRSF4", "TNFRSF18", "LAG3", "VSIR", "ITGB2", "CD2")

p <- DotPlot(TNK, features = genes, cols = c("lightgrey", "#252525"),
             col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6,
             scale.min = 0, scale.max = 50) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("./data/FigS7/Tcells_Ligands_Receptors1.pdf", p, width = 10, height = 10, dpi = 300)

# T cells Ligands and Receptors 2: Complement pathway ligands
ComPClassic <- c("C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C4A", "C4B")
ComPLectin <- c("MBL2", "FCN1", "FCN2", "FCN3", "MASP1", "MASP2")
ComPAlternative <- c("C3", "CFB", "CFD", "CFH", "CFP")
ComPDown <- c("C5", "C6", "C7", "C8A", "C8B", "C9")
ComPLigands <- c(ComPClassic, ComPLectin, ComPAlternative, ComPDown)

p <- DotPlot(TNK, features = ComPLigands, cols = c("lightgrey", "#252525"),
             col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6,
             scale.min = 0, scale.max = 50) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("./data/FigS7/Tcells_Ligands_Receptors2.pdf", p, width = 10, height = 10, dpi = 300)

print("All T/NK cell plots generated successfully!")

# ------------------------------------------------------------------------------
# B Cell Analysis
# ------------------------------------------------------------------------------
# Load B cell data
Bsub <- readRDS("./data/FigS7/S7_Bsub.RDS")

# Set color palette
mypalette <- c("#A6CEE3", "#FB9A99", "#1B9E77", "#1F78B4", "#E31A1C", "#CAB2D6",
               "darkred", "#33A02C", "#FDBF6F", "darkgrey", "#7570B3", "#D95F02",
               "#6A3D9A", "#B15928", "#FF7F00", "#FFFF99", "black", "navy",
               "firebrick", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
               "#FFD92F", "#E5C494")

# B cells Ligands and Receptors 1: Co-stimulatory and co-inhibitory molecules
genes <- c("ICOSLG", "CD72", "LGALS9", "CD40", "CD86", "CD70", "CD58",
           "TNFRSF14", "ICAM1")

p <- DotPlot(Bsub, features = genes, cols = c("lightgrey", "#252525"),
             col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6,
             scale.min = 0, scale.max = 50) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("./data/FigS7/Bcells_Ligands_Receptors1.pdf", p, width = 10, height = 10, dpi = 300)

# B cells Ligands and Receptors 2: Complement pathway ligands
ComPClassic <- c("C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C4A", "C4B")
ComPLectin <- c("MBL2", "FCN1", "FCN2", "FCN3", "MASP1", "MASP2")
ComPAlternative <- c("C3", "CFB", "CFD", "CFH", "CFP")
ComPDown <- c("C5", "C6", "C7", "C8A", "C8B", "C9")
ComPLigands <- c(ComPClassic, ComPLectin, ComPAlternative, ComPDown)

p <- DotPlot(Bsub, features = ComPLigands, cols = c("lightgrey", "#252525"),
             col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6,
             scale.min = 0, scale.max = 50) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("./data/FigS7/Bcells_Ligands_Receptors2.pdf", p, width = 10, height = 10, dpi = 300)

print("All B cell plots generated successfully!")

# ==============================================================================
# End of Script
# ==============================================================================
