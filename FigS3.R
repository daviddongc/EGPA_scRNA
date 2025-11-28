# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S3
# ==============================================================================
# Description: This script generates Supplementary Figure S3
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
library(dplyr)
library(viridis)
library(patchwork)

# Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
granu <- readRDS("./data/FigS3/granu.RDS")
granu.allmarkers <- readRDS("./data/FigS3/granu.allmarkers.rds")

# Define color palette
mypalette <- c("#A6CEE3", "#FB9A99", "#1F78B4", "#E31A1C", "#CAB2D6", "darkred",
               "#33A02C", "#FDBF6F", "darkgrey", "#7570B3", "#D95F02", "#6A3D9A",
               "#B15928", "#1B9E77", "#FF7F00", "#FFFF99", "black", "navy", "firebrick")

# ==============================================================================
# FIGURE S3A: Heatmap of Top Markers
# ==============================================================================
DefaultAssay(granu) <- "SCT"
a <- granu.allmarkers[which(granu.allmarkers[, 2] >= 0.25 & granu.allmarkers[, 5] <= 0.05), ]
top15 <- a %>%
  group_by(cluster) %>%
  slice(1:15)

p <- DoHeatmap(granu, features = top15$gene, disp.min = -2, disp.max = 2)
bks <- seq(-2.1, 2.1, by = 0.1)
heatmapColors <- colorRampPalette(rev(c(brewer.pal(7, "RdBu"))))(length(bks) - 1)
p + scale_fill_gradientn(colors = heatmapColors)

pdf(file = "./data/FigS3/DoHeatmap.pdf", width = 12, height = 10)
print(p + scale_fill_gradientn(colors = heatmapColors))
dev.off()

# ==============================================================================
# FIGURE S3B: SCENIC Transcription Factor Heatmap
# ==============================================================================
AucMTX <- read.csv("./data/FigS3/auc_mtxCopy.csv", sep = "\t", header = TRUE, row.names = 1)
AucMTX <- t(AucMTX)
cells <- unname(unlist(x = CellsByIdentities(granu)))
AucMTX <- AucMTX[, which(colnames(AucMTX) %in% cells)]
sample.selected <- intersect(cells, colnames(AucMTX))
pdata <- granu@meta.data[sample.selected, c("subset_2", "group")]
a <- granu.allmarkers[which(granu.allmarkers[, 2] >= 0.25 & granu.allmarkers[, 5] <= 0.05), ]
a <- a %>% group_by(cluster)
gene.selected <- intersect(a$gene, rownames(AucMTX))
gene.selected <- c("CEBPE", "ZBTB7A", "CUX1", "RXRA", "TCF7L2", "POU2F2", "YBX1",
                   "MAFB", "KLF4", "KLF10", "IRF8", "ATF3", "FOSB", "ETS2", "CREM",
                   "NFE2L2", "FOSL2", "REL", "KLF3", "FOS", "CEBPD", "JUNB", "KLF2",
                   "IRF1", "SPI1", "IRF7", "CEBPB", "ELF1", "BCL6", "STAT1", "TFEC",
                   "MAFF", "GATA1", "ENO1", "ATF4", "KLF5", "GATA2", "AHR", "PBX1",
                   "MITF", "BHLHE40", "EGR1", "JUN")

data <- AucMTX[gene.selected, sample.selected]
data <- t(scale(t(data), scale = TRUE, center = TRUE))
data <- MinMax(data = data, min = -2, max = 2)

pheatmap(data, annotation = pdata, show_colnames = FALSE, cluster_rows = FALSE,
         cluster_cols = FALSE, color = colorRampPalette(inferno(9))(100))

# ==============================================================================
# FIGURE S3C: Violin Plots of Granule Proteins
# ==============================================================================
DefaultAssay(granu) <- "SCT"
specific_granule_genes <- c("MMP9", "HP", "CRISP3", "OLFM4", "LTF", "PGLYRP1",
                            "CHI3L1", "PLBD1")
azurophilic_granule_genes <- c("BPI", "AZU1", "PRTN3", "RNASE3", "CTSG", "MPO",
                               "ELANE", "GUSB", "HEXA")

granu1 <- granu
granu1@active.ident <- factor(granu1@active.ident,
                              levels = rev(levels(granu1@active.ident)))

# Plot granule genes
p111 <- VlnPlot(granu1, slot = "data", assay = "SCT", stack = TRUE,
                features = c(specific_granule_genes, azurophilic_granule_genes),
                flip = FALSE, same.y.lims = FALSE, fill.by = 'ident',
                cols = rev(c(mypalette)[1:14]), pt.size = 0)

p111$layers[[1]]$aes_params$colour <- NA
p111$facet$params$switch <- "x"

p111 <- p111 +
  ggtitle("granu") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(hjust = 0, angle = 270, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    panel.border = element_blank()
  )

print(p111)

# Plot secretory vesicle genes
gene <- c("SERPINA1", "VCAN", "CFP", "ASAH1", "NPC2", "LAMP1", "LAMP2", "GNS",
          "CTSZ", "CTSC", "PRCP", "HEXB", "CTSH")

p112 <- VlnPlot(granu1, slot = "data", assay = "SCT", stack = TRUE,
                features = gene, flip = FALSE, same.y.lims = FALSE,
                fill.by = 'ident', cols = rev(c(mypalette)[1:14]), pt.size = 0)

p112$layers[[1]]$aes_params$colour <- NA
p112$facet$params$switch <- "x"

p112 <- p112 +
  ggtitle("granu") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(hjust = 0, angle = 270, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    panel.border = element_blank()
  )

print(p112)
p111 / p112

# ==============================================================================
# FIGURE S3D: Dotplots of Ligands and Receptors
# ==============================================================================

# CSF ligands and receptors
CSF <- c("CSF1")
CSFR <- c("CSF1R", "CSF2RA", "CSF2RB", "CSF3R")

p1 <- DotPlot(granu, features = rev(c(CSF)), assay = "SCT",
              cols = c("lightgrey", "#252525"), col.min = -1, col.max = 1,
              dot.min = 0, dot.scale = 6, scale.min = 0, scale.max = 20) +
  coord_flip() +
  xlab("Ligands") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

p2 <- DotPlot(granu, features = rev(c(CSFR)), assay = "SCT",
              cols = c("lightgrey", "#252525"), col.min = -1, col.max = 1,
              dot.min = 0, dot.scale = 6, scale.min = 0, scale.max = 20) +
  coord_flip() +
  xlab("Receptors") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

p1 / p2

# Chemokine ligands and receptors
ligands <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL20", "CCL22",
             "CCL23", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL10", "CXCL16", "CXCL17")
receptors <- c("CCR1", "CCR2", "CCR3", "CCR5", "CXCR1", "CXCR2", "CXCR4", "CX3CR1")

p1 <- DotPlot(granu, features = rev(c(ligands)), assay = "SCT",
              cols = c("lightgrey", "#252525"), col.min = -1, col.max = 1,
              dot.min = 0, dot.scale = 6, scale.min = 0, scale.max = 20) +
  coord_flip() +
  xlab("Ligands") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

p2 <- DotPlot(granu, features = rev(c(receptors)), assay = "SCT",
              cols = c("lightgrey", "#252525"), col.min = -1, col.max = 1,
              dot.min = 0, dot.scale = 6, scale.min = 0, scale.max = 20) +
  coord_flip() +
  xlab("Receptors") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

p1 / p2

# Interleukin ligands and receptors
IL_ligands <- c("IL1A", "IL1B", "IL10", "IL15", "IL16", "IL18", "IL23A", "IL32")
IL_receptors <- c("IL1R1", "IL1R2", "IL1RL1", "IL1RN", "IL2RA", "IL2RB", "IL2RG",
                  "IL3RA", "IL4R", "IL5RA", "IL6R", "IL9R", "IL10RA", "IL11RA",
                  "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1",
                  "IL21R", "IL27RA")

p1 <- DotPlot(granu, features = rev(c(IL_ligands)), assay = "SCT",
              cols = c("lightgrey", "#252525"), col.min = -1, col.max = 1,
              dot.min = 0, dot.scale = 6, scale.min = 0, scale.max = 20) +
  coord_flip() +
  xlab("Ligands") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

p2 <- DotPlot(granu, features = rev(c(IL_receptors)), assay = "SCT",
              cols = c("lightgrey", "#252525"), col.min = -1, col.max = 1,
              dot.min = 0, dot.scale = 6, scale.min = 0, scale.max = 20) +
  coord_flip() +
  xlab("Receptors") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

p1 / p2

# ==============================================================================
# End of Script
# ==============================================================================
