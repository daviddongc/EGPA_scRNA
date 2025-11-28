# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S10
# ==============================================================================
# Description: This script generates Supplementary Figure S10
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic
# granulomatosis with polyangiitis> publication
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)
library(sctransform)
library(RColorBrewer)
library(MAST)
library(data.table)
library(rsvd)
library(Rmagic)
library(reticulate)
library(ggsignif)

# Clear workspace
rm(list = ls())

# ==============================================================================
# FIGURE S10A: Cell Type Proportions Analysis
# ==============================================================================

# Load filtered Seurat object
scRNA_harmony_group_filtered <- readRDS("./data/FigS10/F4_i.rds")

# Define custom color palette for cell types
custom_palette <- c(
  "Basal" = "#A6CEE3",
  "Mature ciliated" = "#8B0000",
  "Proximal ciliated" = "#33A02C",
  "Deuterosome-stage ciliated" = "#FDBF6F",
  "Club cells" = "#2178AE",
  "Cycling" = "#FF7F00",
  "Differentiating basal" = "#F49999",
  "Goblet-1" = "#E21B1B",
  "Goblet-2" = "#CAB2D6",
  "Ionocyte" = "#FF00FF",
  "Neuroendocrine" = "#A9A9A9"
)

# Set parameters (should match the original analysis)
dims <- 10
k_param <- 20
resolution <- 1

# Calculate cell proportions by group
cell_proportions <- scRNA_harmony_group_filtered@meta.data %>%
  group_by(group, cell_type) %>%
  tally() %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Ensure that 'group' is a factor and levels are ordered
cell_proportions$group <- factor(cell_proportions$group, levels = c("E", "Remission", "Relapsing"))

# Plot the stacked bar plot with the custom palette
p3 <- ggplot(cell_proportions, aes(x = factor(group), y = proportion, fill = factor(cell_type))) +
  geom_bar(stat = "identity") +
  labs(title = "Cell Type Proportions across Different Groups",
       x = "Group",
       y = "Proportion of Cells",
       fill = "Seurat Cluster") +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(title = "Seurat Clusters"))

# Display plot
print(p3)

# Save plots
ggsave(filename = paste0('./data/FigS10/03.cell_Proportions_dims', dims, '_k', k_param, '_res', resolution, '.png'),
       p3, width = 4, height = 5, dpi = 80)
ggsave(filename = paste0('./data/FigS10/03.cell_Proportions_dims', dims, '_k', k_param, '_res', resolution, '.pdf'),
       p3, width = 4, height = 5, dpi = 300)

cat("Successfully generated 03.cell_Proportions plots\n")

# ==============================================================================
# FIGURE S10C & D: Fibroblast UMAP and Feature Plots
# ==============================================================================

# Load fibro data
fibro <- readRDS("./data/FigS10/S10_fibro.RDS")

# Define color palette
mypalette <- c("#A6CEE3", "#FB9A99", "#1F78B4", "#E31A1C", "#CAB2D6", "darkred",
               "#33A02C", "#FDBF6F", "darkgrey", "#7570B3", "#D95F02", "#6A3D9A",
               "#B15928", "#1B9E77", "#FF7F00", "#FFFF99", "black", "navy", "firebrick")

# Set identities based on subset_2 (cell type annotations)
Idents(fibro) <- factor(fibro$subset_2, levels = c("alvF", "advF", "myofibroblast", "asm", "vsm",
                                                    "pericyte", "schwann", "art", "vein", "bro", "lym"))

# Generate UMAP plot
p <- DimPlot(fibro,
             reduction = "umap",
             cols = mypalette,
             pt.size = 2,
             label = FALSE)

# Save plot
ggsave("./data/FigS10/S10_c.pdf", p, width = 10, height = 10, dpi = 300)

# IGF1 feature plot
p <- FeaturePlot(fibro, features = c("IGF1"),
                 cols = rev(brewer.pal(11, "RdYlBu")[1:7]),
                 reduction = "umap", min.cutoff = 0, max.cutoff = 6,
                 pt.size = 3, combine = T, order = T)
AugmentPlot(p)
ggsave("./data/FigS10/S10_d_IGF1.pdf", p, width = 10, height = 10, dpi = 300)

# CCL11 feature plot
p <- FeaturePlot(fibro, features = c("CCL11"),
                 cols = rev(brewer.pal(11, "RdYlBu")[1:7]),
                 reduction = "umap", min.cutoff = 0, max.cutoff = 6,
                 pt.size = 3, combine = T, order = T)
AugmentPlot(p)
ggsave("./data/FigS10/S10_d_CCL11.pdf", p, width = 10, height = 10, dpi = 300)

# ==============================================================================
# FIGURE S10E: Differential Expression Gene Analysis (advF)
# ==============================================================================

# Load data
fibro <- readRDS("./data/FigS10/S10_fibro.RDS")

# Subset advF cells
spe <- subset(fibro, idents = c("advF"))

# Remove ribosomal and mitochondrial genes
RPS <- grep(pattern = "^RPS", x = row.names(spe[["RNA"]]@data), value = T)
RPL <- grep(pattern = "^RPL", x = row.names(spe[["RNA"]]@data), value = T)
MT <- grep(pattern = "^MT-", x = row.names(spe[["RNA"]]@data), value = T)
spe[["SCT"]]@data <- spe[["SCT"]]@data[-which(row.names(spe[["SCT"]]@data) %in% c(RPS, RPL, MT)), ]

# Convert to MAST object
spe <- FromMatrix(as.matrix(spe[["SCT"]]@data), spe@meta.data)

# Calculate CDR (Cellular Detection Rate)
cdr2 <- colSums(assay(spe) > 0)
qplot(x = cdr2, y = colData(spe)$nFeature_RNA) + xlab('New CDR') + ylab('Old CDR')
colData(spe)$cngeneson <- scale(cdr2)

# Set reference level
colData(spe)$group <- factor(colData(spe)$group, levels = c("C", "A", "E"))

# Note: Need to run zlm analysis first (commented out), then load the saved results
# zlmCond <- zlm(~group, spe)
# summaryCondA <- summary(zlmCond, doLRT = "groupA")
# summaryCondE <- summary(zlmCond, doLRT = "groupE")

# Load previously saved MAST analysis results
summaryCondA <- readRDS("./data/FigS10/S10_advF_Ashma_summaryCond.rds")
summaryCondE <- readRDS("./data/FigS10/S10_advF_EGPA_summaryCond.rds")

# Extract DEGs for Asthma group
summaryDtA <- summaryCondA$datatable
fcHurdleA <- merge(summaryDtA[contrast == "groupA" & component == 'H', .(primerid, `Pr(>Chisq)`)],
                   summaryDtA[contrast == "groupA" & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid')
fcHurdleSigA <- merge(fcHurdleA[`Pr(>Chisq)` <= .05 & abs(coef) >= 0.5], as.data.table(mcols(spe)), by = 'primerid')
setorder(fcHurdleSigA, coef)

# Extract DEGs for EGPA group
summaryDtE <- summaryCondE$datatable
fcHurdleE <- merge(summaryDtE[contrast == "groupE" & component == 'H', .(primerid, `Pr(>Chisq)`)],
                   summaryDtE[contrast == "groupE" & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid')
fcHurdleSigE <- merge(fcHurdleE[`Pr(>Chisq)` <= .05 & abs(coef) >= 0.5], as.data.table(mcols(spe)), by = 'primerid')
setorder(fcHurdleSigE, coef)

# Function for gene labeling
LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0.1, adj.y.t = 0.1, adj.x.s = 0.1,
    adj.y.s = 0.1, text.size = 3, segment.size = 0.1) {
    for (i in genes) {
        x1 <- exp.mat[i, "coef"]
        y1 <- exp.mat[i, "coef.1"]
        plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t,
            label = i, size = text.size)
        plot <- plot + annotate("segment", x = x1, xend = x1 + adj.x.s, y = y1,
            yend = y1 + adj.y.s, size = segment.size)
    }
    return(plot)
}

# Prepare plotting data
rownames(fcHurdleA) <- fcHurdleA$primerid
rownames(fcHurdleE) <- fcHurdleE$primerid
gene <- intersect(fcHurdleA$primerid, fcHurdleE$primerid)
data <- data.frame(cbind(fcHurdleA[gene, ], fcHurdleE[gene, ]))
rownames(data) <- gene

# Filter and color coding
data <- data[which(!is.na(data[, "coef"] & data[, "coef.1"])), ]
data$color <- rep("grey", length(data[, 1]))
data[which(rownames(data) %in% fcHurdleSigA$primerid & rownames(data) %in% fcHurdleSigE$primerid), "color"] <- "both"
data[which(rownames(data) %in% fcHurdleSigA$primerid & !rownames(data) %in% fcHurdleSigE$primerid), "color"] <- "A"
data[which(!rownames(data) %in% fcHurdleSigA$primerid & rownames(data) %in% fcHurdleSigE$primerid), "color"] <- "E"

# Plot fold change scatter plot
sigGenes <- unique(c(fcHurdleSigA$primerid, fcHurdleSigE$primerid))
p <- ggplot(data, aes(x = coef, y = coef.1)) + geom_point(aes(color = color)) + theme_bw()
p <- LabelPoint(p, genes = sigGenes, data)

# Save plot
ggsave("./data/FigS10/advF_DEGs.pdf", p, dpi = 300)

cat("Successfully generated advF_DEGs.pdf\n")

# ==============================================================================
# FIGURE S10F: Adventitial Fibroblast Violin Plot
# ==============================================================================

fibro <- readRDS("./data/FigS10/S10_fibro.RDS")
Adventitial_fibroblast <- subset(fibro, idents = c('advF'))

# NOTE: Configure your own conda environment and Python path before running MAGIC
# conda_path <- "/path/to/your/conda"
# env_path <- "/path/to/your/magic/env"
# use_condaenv(env_path, conda = conda_path, required = TRUE)
# use_python("/path/to/your/python", required = TRUE)
# py_config()

set.seed(123)
expression_matrix <- as.matrix(Adventitial_fibroblast[["SCT"]]@data)
MAGIC_data <- magic(expression_matrix, seed = 123)

Adventitial_fibroblast[["MAGIC_STC"]] <- CreateAssayObject(MAGIC_data$result)

# Set group labels
Adventitial_fibroblast$group <- factor(
  Adventitial_fibroblast$group,
  levels = c("C", "A", "E"),
  labels = c("Control", "Asthma", "EGPA")
)

# Set active assay and ident
Idents(Adventitial_fibroblast) <- "group"
Adventitial_fibroblast@active.assay <- "MAGIC_STC"

# Define comparisons
comparisons <- list(
  c("Control", "Asthma"),
  c("Control", "EGPA"),
  c("Asthma", "EGPA")
)

# Create stacked violin plot
p3 <- VlnPlot(Adventitial_fibroblast,
            stack = TRUE,
            features = c("IGF1", "CCL11"),
            flip = T,
            same.y.lims = FALSE,
            fill.by = 'ident',
            cols = c("#5BBB93", "#377EB8", "#BA2833")) +
  scale_y_continuous(limits = c(0, 0.8)) +
  geom_signif(
    comparisons = comparisons,
    map_signif_level = FALSE,
    textsize = 3,
    test = wilcox.test,
    step_increase = 0.14,
    size = 1,
    tip_length = 0.03,
    family = "Arial",
    y_position = 0.55) +
  ggtitle("Adventitial fibroblast") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())

# Display plot
print(p3)

# Save plot
ggsave(filename = "./data/FigS10/S10_f.pdf", p3, width = 4, height = 6)
ggsave(filename = "./data/FigS10/S10_f.png", p3, width = 4, height = 6, dpi = 600)

# ==============================================================================
# End of Script
# ==============================================================================
