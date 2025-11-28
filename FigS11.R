#!/usr/bin/env Rscript
# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S11
# ==============================================================================
# Description: This script generates Supplementary Figure S11
# IGF1R, IL17RA, IL17RB expression plots
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic
# granulomatosis with polyangiitis> publication
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)

# Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# Setup Output Directory
# ------------------------------------------------------------------------------
output_dir <- "./data/FigS11"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# ------------------------------------------------------------------------------
# Function Definitions
# ------------------------------------------------------------------------------
# Function to save feature plot with specific parameters
save_feature_plot <- function(seurat_obj, gene, output_path, pt.size = 0.5, max.cutoff = NULL, alpha = NULL) {
  cat("Generating plot for", gene, "...\n")

  # Set default assay to SCT (use data slot)
  DefaultAssay(seurat_obj) <- "SCT"

  # Create FeaturePlot with default colors
  p <- FeaturePlot(
    seurat_obj,
    features = gene,
    reduction = "umap",
    slot = "data",
    pt.size = pt.size,
    max.cutoff = max.cutoff
  )

  # Apply alpha transparency if specified
  if (!is.null(alpha)) {
    p$layers[[1]]$aes_params$alpha <- alpha
  }

  # Save as PDF
  ggsave(
    filename = output_path,
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )

  cat("Saved:", output_path, "\n")
}

# ==============================================================================
# FIGURE S11A: Combined Dataset Plots (All Cells)
# ==============================================================================
combined <- readRDS("./data/FigS11/combined.RDS")
cat("Loaded combined.RDS\n")

# IGF1R: max.cutoff = 0.5, alpha = 0.2
save_feature_plot(combined, "IGF1R",
                  file.path(output_dir, "IGF1R_scatter.pdf"),
                  pt.size = 1,
                  max.cutoff = 0.5,
                  alpha = 0.2)

# IL17RA: max.cutoff = "q95"
save_feature_plot(combined, "IL17RA",
                  file.path(output_dir, "IL17RA_scatter.pdf"),
                  pt.size = 1,
                  max.cutoff = "q95")

# IL17RB: max.cutoff = 3
save_feature_plot(combined, "IL17RB",
                  file.path(output_dir, "IL17RB_scatter.pdf"),
                  pt.size = 1,
                  max.cutoff = 3)

rm(combined)
gc()

# ==============================================================================
# FIGURE S11B: Granulocyte Subset Plots
# ==============================================================================
granu <- readRDS("./data/FigS11/granu.RDS")
cat("Loaded granu.RDS\n")

# IGF1R: pt.size = 1, max.cutoff = 0.4
save_feature_plot(granu, "IGF1R",
                  file.path(output_dir, "IGF1R_scatter_granu.pdf"),
                  pt.size = 1,
                  max.cutoff = 0.4)

# IL17RA: pt.size = 1, max.cutoff = "q95"
save_feature_plot(granu, "IL17RA",
                  file.path(output_dir, "IL17RA_granu_scatter.pdf"),
                  pt.size = 1,
                  max.cutoff = "q95")

rm(granu)
gc()

# ==============================================================================
# FIGURE S11C: T & NK Cell Subset Plot
# ==============================================================================
TNK <- readRDS("./data/FigS11/TNK.RDS")
cat("Loaded TNK.RDS\n")

# IL17RB: pt.size = 0.5, max.cutoff = "q95"
save_feature_plot(TNK, "IL17RB",
                  file.path(output_dir, "IL17RB_scatter_TNK.pdf"),
                  pt.size = 0.5,
                  max.cutoff = "q95")

rm(TNK)
gc()

cat("\nAll plots completed successfully!\n")

# ==============================================================================
# End of Script
# ==============================================================================
