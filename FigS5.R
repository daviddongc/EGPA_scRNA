#!/usr/bin/env Rscript
# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S5
# ==============================================================================
# Description: This script generates Supplementary Figure S5
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic
# granulomatosis with polyangiitis> publication
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(reticulate)
library(writexl)
library(ggplot2)
library(readxl)

# Clear workspace
rm(list = ls())

# ==============================================================================
# FIGURE S5E (Left): Differential Expression Heatmap
# ==============================================================================
# Load granulocyte data
granu <- readRDS("./data/FigS5/granu.RDS")
cat("Dataset loaded successfully!\n")

# Comparing 2 subsets and plotting the heatmap
DEGs <- FindMarkers(granu, ident.1 = "cirEos", ident.2 = "Eos", logfc.threshold = 0.1)

# Remove ribosomal and mitochondrial genes
RPgenes <- grep(pattern = "^RP", x = rownames(DEGs), v = T)
MTgenes <- grep(pattern = "^MT-", x = rownames(DEGs), v = T)
DEGs <- DEGs[which(!rownames(DEGs) %in% c(RPgenes, MTgenes)), ]

# Save DEG results
output_dir <- "./data/FigS5"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_xlsx(DEGs, path = file.path(output_dir, "cirEos_vs_Eos_DEGs.xlsx"))

# Select top differential genes
Mast2Genes <- rownames(head(DEGs[order(DEGs[, 2]), ], 50))      # top 50 up in Eos
Mast1Genes <- rev(rownames(tail(DEGs[order(DEGs[, 2]), ], 63))) # top 63 up in cirEos
genes <- c(Mast1Genes, Mast2Genes)

# Extract cells
celltag1 <- rownames(granu@meta.data[which(granu@meta.data[, "subset_2"] %in% c("cirEos")), ])
celltag2 <- rownames(granu@meta.data[which(granu@meta.data[, "subset_2"] %in% c("Eos")), ])

# Downsample using geosketch for visualization
pc <- Embeddings(granu, reduction = "pca")
geosketch <- import("geosketch")
sketch.fraction <- 0.3
sketch.size <- as.integer(dim(granu)[2] * sketch.fraction)
sketch.indice <- geosketch$gs(pc, sketch.size)
sample.selected <- colnames(granu)[unlist(sketch.indice)]
sample.selected <- intersect(c(celltag1, celltag2), sample.selected)

# Remove interfering genes (HLA, interferon-induced, etc.)
genes <- genes[which(!genes %in% c(
  "HLA-B", "IFI6", "HLA-A", "HLA-DRA", "RSAD2", "HLA-DRB1",
  "HLA-DPA1", "CTSB", "B2M", "HLA-DPB1", "CD74", "SLPI", "LYZ"
))]

# Extract and normalize data
data <- as.matrix(granu[["SCT"]]@data[unique(genes), sample.selected])
data <- t(scale(t(data), scale = TRUE, center = TRUE))

# Function to scale values between min and max
MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

data <- MinMax(data = data, min = -2, max = 2)

# Plot heatmap
pdf(file.path(output_dir, "cirEos_vs_Eos_heatmap.pdf"), width = 8, height = 12)
pheatmap(data,
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(rev(c(brewer.pal(9, "RdBu"))))(100)
)
dev.off()

cat("Heatmap saved to:", file.path(output_dir, "cirEos_vs_Eos_heatmap.pdf"), "\n")

# ==============================================================================
# FIGURE S5E (Right): GO Enrichment Bar Plots
# ==============================================================================

# Plotting the Metascape GO barplots
# Read GO enrichment results from Excel (Sheet 2 contains GO terms)
# This file should have GO results from Metascape
GO <- read_excel(file.path(output_dir, "cirEos_vs_Eos_DEGs.xlsx"), sheet = 2)

# Split GO results: first 8 rows for cirEos, rows 9-16 for Eos
GO1 <- GO[1:8, ]   # cirEos upregulated genes' GO terms
GO2 <- GO[9:16, ]  # Eos upregulated genes' GO terms

# Define colors
color1 <- "#7570AF"  # Purple for cirEos
color2 <- "#D95F0D"  # Orange for Eos

# Plot cirEos GO bar plot
p1 <- ggplot(GO1, aes(x = x, y = pvalue, fill = pvalue)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = color1) +
  labs(
    x = "GO Terms",
    y = "-log10(p-value)",
    title = "Circulatory Eosinophils"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot Eos GO bar plot
p2 <- ggplot(GO2, aes(x = x, y = pvalue, fill = pvalue)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_gradient(low = "lightgrey", high = color2) +
  labs(
    x = "GO Terms",
    y = "-log10(p-value)",
    title = "Tissue Eosinophils"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plots
ggsave(file.path(output_dir, "cirEos_bar.pdf"),
  plot = p1, width = 8, height = 6, dpi = 300
)
ggsave(file.path(output_dir, "Eos_bar.pdf"),
  plot = p2, width = 8, height = 6, dpi = 300
)

cat("GO bar plots saved to:\n")
cat("  -", file.path(output_dir, "cirEos_bar.pdf"), "\n")
cat("  -", file.path(output_dir, "Eos_bar.pdf"), "\n")

# ==============================================================================
# End of Script
# ==============================================================================

