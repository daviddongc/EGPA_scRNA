# ==============================================================================
# EGPA scRNA-seq Analysis: Figure 3
# ==============================================================================
# Description: This script generates Figure 3A, 3B, 3C, 3F
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic
# granulomatosis with polyangiitis> publication
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
library(velocyto.R)
library(stringr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reticulate)

# Clear workspace
rm(list = ls())

# ==============================================================================
# FIGURE 3A: RNA Velocity Analysis
# ==============================================================================

# Load metadata
Bsub_meta <- readRDS("./data/Fig3/Bsub_meta.rds")

# Get list of loom files
listdir <- dir("./data/Fig3/loom_72")

# Step 1: Generating the combined ematMat / nmatMat / ambMat
emat <- list()
nmat <- list()
amb <- list()

for (i in 1:length(listdir)) {
  data <- read.loom.matrices(paste("./data/Fig3/loom_72", listdir[i], sep = "/"))
  emat[i] <- data$spliced
  nmat[i] <- data$unspliced
  amb[i] <- data$ambiguous
}

ematMat <- data.frame()
nmatMat <- data.frame()
ambMat <- data.frame()

for (i in 1:length(listdir)) {
  ldat <- emat[[i]]
  colnames <- strsplit(colnames(ldat), split = ":")
  frame <- matrix(unlist(colnames), ncol = 2, byrow = TRUE)
  frame[, 2] <- unlist(strsplit(frame[, 2], split = "x"))
  colnames <- paste(frame[, 1], frame[, 2], sep = "_")

  colnames(emat[[i]]) <- colnames
  colnames(nmat[[i]]) <- colnames
  colnames(amb[[i]]) <- colnames

  inter <- intersect(rownames(Bsub_meta), colnames)
  if (i == 1) {
    ematMat <- emat[[i]][, inter]
    nmatMat <- nmat[[i]][, inter]
    ambMat <- amb[[i]][, inter]
  } else {
    ematMat <- cbind(ematMat, emat[[i]][, inter])
    nmatMat <- cbind(nmatMat, nmat[[i]][, inter])
    ambMat <- cbind(ambMat, amb[[i]][, inter])
  }
}

# Save matrices
saveRDS(ematMat, file = "./data/Fig3/ematMat.rds")
saveRDS(nmatMat, file = "./data/Fig3/nmatMat.rds")
saveRDS(ambMat, file = "./data/Fig3/ambMat.rds")

# Step 2: Plotting the velocity
ematMat <- readRDS("./data/Fig3/ematMat.rds")
nmatMat <- readRDS("./data/Fig3/nmatMat.rds")

# Load B cell subset data
Bsub <- readRDS("./data/Fig3/Bsub.RDS")
Bsub <- subset(Bsub, cells = colnames(ematMat))

# Perform geometric sketching
pc <- Embeddings(Bsub, reduction = "pca")
geosketch <- import("geosketch", convert = TRUE)
sketch.fraction <- 0.5
sketch.size <- as.integer(dim(Bsub)[2] * sketch.fraction)
sketch.indice <- geosketch$gs(pc, sketch.size)
cell.selected <- colnames(Bsub)[unlist(sketch.indice)]

Bsub <- subset(Bsub, cells = cell.selected)
cluster <- Bsub$subset_2
emb <- Embeddings(Bsub, reduction = "umap")
ematMat <- ematMat[, cell.selected]
nmatMat <- nmatMat[, cell.selected]

# Create UMAP plot
mypalette <- c("#A6CEE3", "#FB9A99", "#1B9E77", "#1F78B4")
pdf("./data/Fig3/Bsub_UMAP.pdf")
DimPlot(Bsub, reduction = "umap", cols = mypalette, pt.size = 1, label = FALSE)
dev.off()

# Calculate velocity
cell.dist <- as.dist(1 - armaCor(t(Embeddings(Bsub, reduction = "pca"))))
ematMat <- filter.genes.by.cluster.expression(ematMat, cluster, min.max.cluster.average = 0.5)
nmatMat <- filter.genes.by.cluster.expression(nmatMat, cluster, min.max.cluster.average = 0.05)
length(intersect(rownames(ematMat), rownames(ematMat)))
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(ematMat, nmatMat, deltaT = 1, kCells = 20,
                                            cell.dist = cell.dist, fit.quantile = fit.quantile)

# Plot velocity
pdf("./data/Fig3/Bsub_velocity.pdf")
show.velocity.on.embedding.cor(emb, rvel.cd, n = 300, scale = 'sqrt', cex = 0.8,
                               arrow.scale = 5, show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                               do.par = FALSE, cell.border.alpha = 0.1, n.cores = 1)
dev.off()

# ==============================================================================
# FIGURE 3B: Marker Gene Expression Feature Plots
# ==============================================================================

# Load B&Plasma cell data
BP <- readRDS("./data/Fig3/Bsub.RDS")
cat("Dataset loaded successfully!\n")
cat("Cells:", ncol(BP), "\n")
cat("Genes:", nrow(BP), "\n")

DefaultAssay(BP) <- "MAGIC_SCT"
pt.size <- 0.5

# Marker genes for B and Plasma cell subsets
features <- c("AICDA", "CD19", "CD27", "CD38", "CR2", "FCRL4", "FCRL5",
              "IGHA1", "IGHD", "IGHG1", "IGHM", "ITGAX", "MS4A1", "SDC1",
              "TNFRSF13B")

colors <- c("grey90", "yellow", "orange", "red", "darkred")

# Create output directory
output_dir <- "./data/Fig3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (gene in features) {
  cat("Plotting", gene, "...\n")

  p <- FeaturePlot(BP,
                   features = gene,
                   reduction = "umap",
                   slot = "data",
                   pt.size = pt.size,
                   cols = colors,
                   max.cutoff = 'q95') +
    ggtitle(gene)

  # Save figure
  output_path <- file.path(output_dir, paste0(gene, ".pdf"))
  ggsave(output_path, plot = p, width = 6, height = 5, dpi = 300)

  cat("Saved to", output_path, "\n")
}

cat("Output directory:", output_dir, "\n")

# ==============================================================================
# FIGURE 3C: Cellular Source Proportions
# ==============================================================================

# Load B&Plasma cell data
Bsub <- readRDS("./data/Fig3/Bsub.RDS")
cat("Dataset loaded successfully!\n")
cat("Cells:", ncol(Bsub), "\n")

# Ensure cell identities are set
Idents(Bsub) <- factor(Bsub$subset_2,
                       levels = c("naive B", "memory B", "FCRL4 memory B", "plasma"))

# Calculate the proportions of different individuals in each subtype
pdata <- data.frame(table(Bsub$subject_id, Bsub$source, Bsub@active.ident))
pdata$TotalNumOfCelltype <- rep(unname(table(Bsub@active.ident)), each = 108)

a <- unique(Bsub@meta.data[, c("subject_id", "group")])
a <- a[order(a$subject_id), ]
pdata$group <- rep(a[, "group"], 16)

pdata$Proportion <- pdata$Freq / pdata$TotalNumOfCelltype
pdata$x <- paste(pdata$Var2, pdata$group, sep = "-")
pdata$x <- factor(pdata$x, levels = c(
  "PBMC-C", "PBMC-A", "PBMC-E",
  "Biopsy-C", "Biopsy-A", "Biopsy-E",
  "BALF-C", "BALF-A", "BALF-E",
  "BALF_EOS-C", "BALF_EOS-A"
))

# Create the stacked bar plot
p <- ggplot(pdata, aes(x = Var3, y = Proportion, fill = x)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()

# Save output
output_path <- file.path(output_dir, "Proportions_of_CellCluster.pdf")
ggsave(output_path, plot = p, width = 10, height = 6, dpi = 300)

cat("Plot saved to:", output_path, "\n")

# ==============================================================================
# FIGURE 3F: Fc Receptor Expression in Granulocytes
# ==============================================================================

# Load granulocyte data
granu <- readRDS("./data/Fig3/granu.RDS")
cat("Dataset loaded successfully!\n")
cat("Cells:", ncol(granu), "\n")

Fcers <- c("FCAR", "TNFSF13B", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A",
           "FCGR3B", "FCGRT", "TRIM21", "FCER1A")

# Plotting the Dotplots of DEGs among different groups
granu$group_cluster <- paste(granu$group, granu$subset_2, sep = "_")
Idents(granu) <- "group_cluster"
granu <- subset(granu, idents = c(
  "C_granu pro1", "A_granu pro1", "E_granu pro1", "C_granu pro2", "A_granu pro2", "E_granu pro2",
  "C_cirNeu-CSF1R", "A_cirNeu-CSF1R", "E_cirNeu-CSF1R", "C_cirNeu-VCN", "A_cirNeu-VCN", "E_cirNeu-VCN",
  "C_cirNeu-TNF", "A_cirNeu-TNF", "E_cirNeu-TNF", "C_Neu-INSIG1", "A_Neu-INSIG1", "E_Neu-INSIG1",
  "C_cirNeu-CXCR2", "A_cirNeu-CXCR2", "E_cirNeu-CXCR2", "C_Neu-IL1RN", "A_Neu-IL1RN", "E_Neu-IL1RN",
  "C_Neu-CSF1", "A_Neu-CSF1", "E_Neu-CSF1", "C_cirEos", "A_cirEos", "E_cirEos",
  "C_Eos", "A_Eos", "E_Eos", "C_Baso", "A_Baso", "E_Baso",
  "C_Mast-ACP5", "A_Mast-ACP5", "E_Mast-ACP5", "C_Mast-CSF1", "A_Mast-CSF1", "E_Mast-CSF1"))

granu@active.ident <- factor(granu@active.ident, levels = rev(c(
  "C_granu pro1", "A_granu pro1", "E_granu pro1", "C_granu pro2", "A_granu pro2", "E_granu pro2",
  "C_cirNeu-CSF1R", "A_cirNeu-CSF1R", "E_cirNeu-CSF1R", "C_cirNeu-VCN", "A_cirNeu-VCN", "E_cirNeu-VCN",
  "C_cirNeu-TNF", "A_cirNeu-TNF", "E_cirNeu-TNF", "C_Neu-INSIG1", "A_Neu-INSIG1", "E_Neu-INSIG1",
  "C_cirNeu-CXCR2", "A_cirNeu-CXCR2", "E_cirNeu-CXCR2", "C_Neu-IL1RN", "A_Neu-IL1RN", "E_Neu-IL1RN",
  "C_Neu-CSF1", "A_Neu-CSF1", "E_Neu-CSF1", "C_cirEos", "A_cirEos", "E_cirEos",
  "C_Eos", "A_Eos", "E_Eos", "C_Baso", "A_Baso", "E_Baso",
  "C_Mast-ACP5", "A_Mast-ACP5", "E_Mast-ACP5", "C_Mast-CSF1", "A_Mast-CSF1", "E_Mast-CSF1")))

P <- DotPlot(granu, features = Fcers, col.min = -1, col.max = 1, dot.min = 0,
             dot.scale = 6, scale.min = 0, scale.max = 100) +
  scale_color_gradientn(colors = colorRampPalette(rev(c(brewer.pal(7, "RdBu"))))(10))

# Save output
output_path <- file.path(output_dir, "Violin_of_CellCluster.pdf")
ggsave(output_path, plot = P, width = 10, height = 6, dpi = 300)

cat("Plot saved to:", output_path, "\n")

# ==============================================================================
# End of Script
# ==============================================================================
