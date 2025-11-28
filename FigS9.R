# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S9
# ==============================================================================
# Description: This script generates Supplementary Figure S9
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic
# granulomatosis with polyangiitis> publication
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(Seurat)
library(Matrix)
library(stringr)
library(harmony)
library(biomaRt)
library(GEOquery)
library(tidyverse)
library(data.table)
library(DoubletFinder)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(ggplot2)
library(ggsci)
library(SeuratExtend)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(patchwork)
library(reshape2)
library(scRNAtoolVis)
library(ggfun)
library(presto)
library(enrichR)
library(ggpubr)
library(magrittr)
library(Rmagic)
library(reticulate)
library(readxl)

# ==============================================================================
# DATA PREPARATION
# ==============================================================================

# ------------------------------------------------------------------------------
# Load and Prepare Macrophage Data
# ------------------------------------------------------------------------------
# Read macrophage data from treatment samples
macroTreat <- readRDS("./data/FigS9/macroTreat.RDS")

# Remove unwanted clusters
macroTreat <- subset(macroTreat, idents = c(18, 17), invert = T)  # Remove VCAN+
macroTreat <- subset(macroTreat, idents = c(6, 30), invert = T)   # Remove MT-+

# Assign cell type annotations
macroTreat$subset_2 <- "macroTreat"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(27))] <- "pDC"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(24))] <- "mDC"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(28, 29))] <- "mitotic macro"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(22))] <- "GPR34 DC"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(21))] <- "EREG DC"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(25))] <- "TREM2 DC"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(0, 11, 12, 23, 19, 9, 15, 7, 10, 20, 16, 26))] <- "macro1"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(3, 4, 14))] <- "macro2"
macroTreat$subset_2[which(Idents(macroTreat) %in% c(1, 5, 8, 2, 13))] <- "IGF1 macro3"

Idents(macroTreat) <- factor(macroTreat$subset_2,
                             levels = c("macro1", "macro2", "IGF1 macro3", "mitotic macro",
                                       "mDC", "GPR34 DC", "EREG DC", "TREM2 DC", "pDC"))
macroTreat$group2 <- "Treat"

# Read baseline macrophage data
macro <- readRDS("./data/FigS9/macro.RDS")
macro$group2 <- "Without"

# Merge baseline and treatment datasets
scRNA <- merge(macro, macroTreat)
system.time(save(scRNA, file = "./data/FigS9/scRNA.orig.rdata"))

# ------------------------------------------------------------------------------
# Process Combined Dataset
# ------------------------------------------------------------------------------
rm(list = ls())
load("./data/FigS9/scRNA.orig.rdata")

# Increase maximum allowed global variable size
options(future.globals.maxSize = 10 * 1024^3)  # Set to 10 GB

DefaultAssay(scRNA) <- "RNA"

# Remove unnecessary assays
scRNA[["SCT"]] <- NULL
scRNA[["MAGIC_SCT"]] <- NULL

# Subset to EGPA, Remission, and Relapsing groups
scRNA_ERR <- subset(scRNA, subset = group %in% c('E', 'Remission', 'Relapsing'))

# Normalization and feature selection
scRNA_ERR <- SCTransform(scRNA_ERR)
scRNA_ERR <- ScaleData(scRNA_ERR)
scRNA_ERR <- FindVariableFeatures(scRNA_ERR)
RunPCA(scRNA_ERR, npcs = 100) -> scRNA_ERR
system.time(save(scRNA_ERR, file = "scRNA_ERR.rdata"))

# ------------------------------------------------------------------------------
# Harmony Batch Correction
# ------------------------------------------------------------------------------

# Harmony correction by group
scRNA_harmony_group <- RunHarmony(scRNA_ERR,
                                  reduction = "pca",
                                  group.by.vars = "group",
                                  reduction.save = "harmony",
                                  num.cores = 55)
FindNeighbors(scRNA_harmony_group, dims = c(1:37), reduction = "harmony") -> scRNA_harmony_group
FindClusters(scRNA_harmony_group, resolution = 1, reduction = "harmony") -> scRNA_harmony_group
RunUMAP(scRNA_harmony_group, dims = c(1:37), reduction = "harmony") -> scRNA_harmony_group
system.time(save(scRNA_harmony_group, file = "scRNA.harmony.group-ERR.rdata"))

# Harmony correction by group2 (treatment status)
scRNA_harmony_group2 <- RunHarmony(scRNA_ERR,
                                   reduction = "pca",
                                   group.by.vars = "group2",
                                   reduction.save = "harmony",
                                   num.cores = 55)
FindNeighbors(scRNA_harmony_group2, dims = c(1:19), reduction = "harmony") -> scRNA_harmony_group2
FindClusters(scRNA_harmony_group2, resolution = 1, reduction = "harmony") -> scRNA_harmony_group2
RunUMAP(scRNA_harmony_group2, dims = c(1:19), reduction = "harmony") -> scRNA_harmony_group2
system.time(save(scRNA_harmony_group2, file = "scRNA.harmony.group2-ERR.rdata"))

# Harmony correction by original identity
scRNA_harmony_orig <- RunHarmony(scRNA_ERR,
                                 reduction = "pca",
                                 group.by.vars = "orig.ident",
                                 reduction.save = "harmony",
                                 num.cores = 55)
FindNeighbors(scRNA_harmony_orig, dims = c(1:37), reduction = "harmony") -> scRNA_harmony_orig
FindClusters(scRNA_harmony_orig, resolution = 1, reduction = "harmony") -> scRNA_harmony_orig
RunUMAP(scRNA_harmony_orig, dims = c(1:37), reduction = "harmony") -> scRNA_harmony_orig
system.time(save(scRNA_harmony_orig, file = "scRNA.harmony.orig-ERR.rdata"))

# ------------------------------------------------------------------------------
# Differential Expression Analysis
# ------------------------------------------------------------------------------

# Find markers using presto
table(scRNA_harmony_group2$SCT_snn_res.1)
system.time(dge <- wilcoxauc(scRNA_harmony_group2, group_by = 'SCT_snn_res.1'))
dge_0.5 <- dge[which(dge$logFC > 0.05 | dge$logFC < -0.5), ]
write.table(dge_0.5, file = './data/FigS9/02.AllMarkers_presto.csv', sep = ',', row.names = F, quote = F)

# Find markers using Seurat
Idents(scRNA_harmony_group2) <- scRNA_harmony_group2$SCT_snn_res.1
system.time(dge1 <- FindAllMarkers(scRNA_harmony_group2, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = 0.25, test.use = 'wilcox',
                                   return.thresh = 0.01))
write.table(dge1, file = './data/FigS9/02.AllMarkers_findmarker.csv', sep = ',', row.names = F, quote = F)

# Remove unwanted clusters
scRNA_harmony_group2_filtered <- subset(scRNA_harmony_group2,
                                        subset = seurat_clusters %in% c(6, 25),
                                        invert = TRUE)

# ------------------------------------------------------------------------------
# Cell Type Annotation
# ------------------------------------------------------------------------------
cluster.ids <- c(
  "0" = "Macrophage1",
  "1" = "IGF1 Macrophage",
  "2" = "Macrophage2",
  "3" = "Macrophage1",
  "4" = "Macrophage1",
  "5" = "Macrophage2",
  "7" = "Macrophage1",
  "8" = "EREG DC",
  "9" = "Macrophage1",
  "10" = "Macrophage1",
  "11" = "GPR34 DC",
  "12" = "mDC2",
  "13" = "Macrophage2",
  "14" = "Macrophage2",
  "15" = "mDC1",
  "16" = "Macrophage1",
  "17" = "GPR34 DC",
  "18" = "Mitotic macrophage",
  "19" = "TREM2 DC",
  "20" = "DC (in PBMC)",
  "21" = "pDC",
  "22" = "Macrophage1",
  "23" = "Macrophage2",
  "24" = "XCXR1 DC"
)

# Assign cell types
scRNA_harmony_group2_filtered <- RenameIdents(scRNA_harmony_group2_filtered, cluster.ids)
scRNA_harmony_group2_filtered$celltype <- Idents(scRNA_harmony_group2_filtered)

# Check cell counts per cell type
table(scRNA_harmony_group2_filtered@meta.data$celltype)

# Set cell type factor levels
scRNA_harmony_group2_filtered$celltype <- factor(
  scRNA_harmony_group2_filtered$celltype,
  levels = c("Macrophage1", "Macrophage2", "IGF1 Macrophage", "Mitotic macrophage",
            "mDC1", "mDC2", "XCXR1 DC", "GPR34 DC", "EREG DC", "TREM2 DC",
            "pDC", "DC (in PBMC)")
)
Idents(scRNA_harmony_group2_filtered) <- "celltype"

# ------------------------------------------------------------------------------
# Update Metadata
# ------------------------------------------------------------------------------
mtd <- scRNA_harmony_group2_filtered@meta.data
mtd <- mtd %>%
  mutate(source = ifelse(!is.na(source) & source != "", source, type))
scRNA_harmony_group2_filtered@meta.data <- mtd

# Recode group labels
scRNA_harmony_group2_filtered$group <- factor(
  scRNA_harmony_group2_filtered$group,
  levels = c("E", "Remission", "Relapsing"),
  labels = c("EGPA", "Cs-remission", "Cs-relapsing")
)

scRNA_harmony_group2_filtered$group2 <- factor(
  scRNA_harmony_group2_filtered$group2,
  levels = c("Without", "Treat"),
  labels = c("Without treatment", "With treatment")
)
system.time(save(scRNA_harmony_group2_filtered, file = "./data/FigS9/scRNA_harmony_group2_filtered.rdata"))
# ==============================================================================
# FIGURE S9B: UMAP of After-Treatment Data
# ==============================================================================

# Load after-treatment data
afterTreat <- readRDS("./data/FigS9/EGPA_afterTreatment.RDS")

# Add metadata
meta <- read.csv("./data/FigS9/metadata.txt", sep = "\t", header = T)
afterTreat$group <- plyr::mapvalues(afterTreat$orig.ident, from = meta$id, to = meta$group)
afterTreat$sampleid <- plyr::mapvalues(afterTreat$orig.ident, from = meta$id, to = meta$sampleid)
afterTreat$type <- plyr::mapvalues(afterTreat$orig.ident, from = meta$id, to = meta$type)

# ------------------------------------------------------------------------------
# Process After-Treatment Data
# ------------------------------------------------------------------------------
afterTreat <- SCTransform(afterTreat)
NormalizeData(afterTreat) -> afterTreat
ScaleData(afterTreat) -> afterTreat
FindVariableFeatures(afterTreat) -> afterTreat
RunPCA(afterTreat, npcs = 50) -> afterTreat
FindNeighbors(afterTreat) -> afterTreat
FindClusters(afterTreat, resolution = 2) -> afterTreat
RunUMAP(afterTreat, dims = 1:15) -> afterTreat

# ------------------------------------------------------------------------------
# Split Cluster 33 (B cells and Mast cells)
# ------------------------------------------------------------------------------
cellname1 <- rownames(afterTreat@meta.data[which(afterTreat@meta.data$seurat_clusters == 33), ])
cellname2 <- names(which(afterTreat[["RNA"]]@data["CD79A", ] > 0))
Bcelltag <- intersect(cellname1, cellname2)

cellname3 <- names(which(afterTreat[["RNA"]]@data["MS4A2", ] > 0))
Mastcelltag <- intersect(cellname1, cellname3)

afterTreat$seurat_clusters <- as.character(afterTreat$seurat_clusters)
afterTreat@meta.data[Bcelltag, "seurat_clusters"] <- 47
afterTreat@meta.data[Mastcelltag, "seurat_clusters"] <- 48
Idents(afterTreat) <- "seurat_clusters"

# Remove EPCAM+ CD68+ double positive cells
cellRe1 <- rownames(afterTreat@meta.data[which(Idents(afterTreat) %in% c(40)), ])
afterTreat <- subset(afterTreat, cells = cellRe1, invert = T)

# ------------------------------------------------------------------------------
# Cell Type Annotation for After-Treatment Data
# ------------------------------------------------------------------------------
afterTreat$coarse_cluster.r1 <- as.vector(Idents(afterTreat))
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(47))] <- "B cells"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(31))] <- "T & NK cells"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in%
                                    c(4, 26, 15, 0, 7, 10, 16, 5, 42, 18, 37, 25, 34, 45))] <- "Macrophage"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in%
                                    c(1, 2, 3, 9, 13, 46, 43, 6, 24, 20, 23, 17, 22, 29, 41, 28, 30, 39, 38, 14))] <- "Epithelial"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(36))] <- "Neutrophils-1"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(12, 44))] <- "Mono DCs"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(48))] <- "Mast cells"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(35))] <- "pDC"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(33))] <- "Eosinophils"
afterTreat$coarse_cluster.r1[which(afterTreat$seurat_clusters %in% c(8, 11, 27, 19, 21, 32))] <- "Neutrophils-2"

Idents(afterTreat) <- "coarse_cluster.r1"
afterTreat@active.ident <- factor(afterTreat@active.ident,
                                  levels = c("Epithelial", "T & NK cells", "B cells",
                                            "Macrophage", "Mono DCs", "pDC", "Neutrophils-1",
                                            "Neutrophils-2", "Eosinophils", "Mast cells"))

# Define color palette
mypalette <- c("#A6CEE3", "#FB9A99", "#1F78B4", "#E31A1C", "#CAB2D6", "darkred",
               "#33A02C", "#FDBF6F", "darkgrey", "#7570B3", "#D95F02", "#6A3D9A",
               "#B15928", "#1B9E77", "#FF7F00", "navy", "black", "firebrick", "#FFFF99")

# Create UMAP plot
p <- DimPlot(afterTreat, reduction = "umap", cols = mypalette, pt.size = 1, label = F)
print(p)

# ==============================================================================
# FIGURE S9D: Violin Plots for Neutrophils-2 and Eosinophils
# ==============================================================================

# ------------------------------------------------------------------------------
# Neutrophils-2 Analysis
# ------------------------------------------------------------------------------
combined_filtered <- readRDS("./data/FigS9/combined_filtered.RDS")
combined_Neutrophils2 <- subset(combined_filtered,
                                subset = coarse_cluster.r1 %in% c("Neutrophils-2"))
combined_Neutrophils2_matrix <- as.matrix(combined_Neutrophils2[["SCT"]]@data)

# NOTE: Configure your own conda environment and Python path before running MAGIC
# conda_path <- "/path/to/your/conda"
# env_path <- "/path/to/your/magic/env"
# use_condaenv(env_path, conda = conda_path, required = TRUE)
# use_python("/path/to/your/python", required = TRUE)
# py_config()

# Run MAGIC imputation
MAGIC_data <- magic(combined_Neutrophils2_matrix, n.jobs = -30)
combined_Neutrophils2[["MAGIC_STC"]] <- CreateAssayObject(MAGIC_data$result)

Idents(combined_Neutrophils2) <- "group"
combined_Neutrophils2@active.assay <- "MAGIC_STC"

plots <- VlnPlot(combined_Neutrophils2,
                 features = c("IRF7", "IFIT1", "ISG15", "MX1", "TNFSF13B",
                             "FCGR1A", "FCGRT", "TRIM21", "IL1B"),
                 pt.size = 0,
                 combine = FALSE,
                 adjust = 1)
CombinePlots(plots = plots, ncol = 1)

p3 <- VlnPlot(combined_Neutrophils2,
              stack = TRUE,
              features = c("IRF7", "IFIT1", "ISG15", "MX1", "TNFSF13B",
                          "FCGR1A", "FCGRT", "TRIM21", "IL1B"),
              flip = T,
              same.y.lims = FALSE,
              fill.by = 'ident',
              cols = c("#BA2833", "#F5A122", "#F19AB8")) +
  scale_y_continuous(limits = c(0, 2)) +
  ggtitle("pDC") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank()
  )
print(p3)

# Process after-treatment Neutrophils-2
afterTreat_Neutrophils2 <- subset(afterTreat_filtered,
                                  subset = coarse_cluster.r1 %in% c("Neutrophils-2"))
afterTreat_Neutrophils2_matrix <- as.matrix(afterTreat_Neutrophils2[["SCT"]]@data)
afterTreat_Neutrophils2_MAGIC_data <- magic(afterTreat_Neutrophils2_matrix, n.jobs = -30)

# Process all Neutrophils-2
all_Neutrophils2 <- subset(scRNA_filtered,
                           subset = coarse_cluster.r1 %in% c("Neutrophils-2"))
all_Neutrophils2_matrix <- as.matrix(all_Neutrophils2[["SCT"]]@data)
all_Neutrophils2_MAGIC_data <- magic(all_Neutrophils2_matrix, n.jobs = -30)
all_Neutrophils2[["MAGIC_STC"]] <- CreateAssayObject(all_Neutrophils2_MAGIC_data$result)

Idents(all_Neutrophils2) <- "group"
all_Neutrophils2@active.assay <- "MAGIC_STC"

p3 <- VlnPlot(all_Neutrophils2,
              stack = TRUE,
              features = c("IRF7", "IFIT1", "ISG15", "MX1", "TNFSF13B",
                          "FCGR1A", "FCGRT", "TRIM21", "IL1B"),
              flip = T,
              same.y.lims = FALSE,
              fill.by = 'ident',
              cols = c("#BA2833", "#F5A122", "#F19AB8")) +
  scale_y_continuous(limits = c(0, 2)) +
  ggtitle("Neutrophils-2") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
print(p3)

ggsave("./data/FigS9/Neutrophils2_VlnPlot.pdf", p3,
       width = 2.5, height = 5, bg = "white")
ggsave("./data/FigS9/Neutrophils2_VlnPlot.png", p3,
       width = 2.5, height = 5, bg = "white", dpi = 600, units = "in")

# ------------------------------------------------------------------------------
# Eosinophils Analysis
# ------------------------------------------------------------------------------
combined_Eosinophils <- subset(combined_filtered,
                               subset = coarse_cluster.r1 %in% c("Eosinophils"))
combined_Eosinophils_matrix <- as.matrix(combined_Eosinophils[["SCT"]]@data)
combined_Eosinophils_MAGIC_data <- magic(combined_Eosinophils_matrix, n.jobs = -30)

afterTreat_Eosinophils <- subset(afterTreat_filtered,
                                 subset = coarse_cluster.r1 %in% c("Eosinophils"))
afterTreat_Eosinophils_matrix <- as.matrix(afterTreat_Eosinophils[["SCT"]]@data)
afterTreat_Eosinophils_MAGIC_data <- magic(afterTreat_Eosinophils_matrix, n.jobs = -30)

all_Eosinophils <- subset(scRNA_filtered,
                          subset = coarse_cluster.r1 %in% c("Eosinophils"))
all_Eosinophils_matrix <- as.matrix(all_Eosinophils[["SCT"]]@data)
all_Eosinophils_MAGIC_data <- magic(all_Eosinophils_matrix, n.jobs = -30)
all_Eosinophils[["MAGIC_STC"]] <- CreateAssayObject(all_Eosinophils_MAGIC_data$result)

Idents(all_Eosinophils) <- "group"
all_Eosinophils@active.assay <- "MAGIC_STC"

p3 <- VlnPlot(all_Eosinophils,
              stack = TRUE,
              features = c("SIGLEC8", "PLAUR", "CCR3", "IL5RA", "IL1RL1",
                          "IL3RA", "IL4R", "IL13RA1"),
              flip = T,
              same.y.lims = FALSE,
              fill.by = 'ident',
              cols = c("#BA2833", "#F5A122", "#F19AB8")) +
  scale_y_continuous(limits = c(0, 1.45)) +
  ggtitle("Eosinophils") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
print(p3)

ggsave("./data/FigS9/Eosinophils_VlnPlot.pdf", p3,
       width = 2.2, height = 4.5, bg = "white")
ggsave("./data/FigS9/Eosinophils_VlnPlot.png", p3,
       width = 2.2, height = 4.5, bg = "white", dpi = 600, units = "in")

# ==============================================================================
# FIGURE S9E: Cell Proportion Bar Plots
# ==============================================================================

# ------------------------------------------------------------------------------
# Proportion by Treatment Status (group2)
# ------------------------------------------------------------------------------
load("./data/FigS9/scRNA_harmony_group2_filtered.rdata")

table_samples_by_clusters <- scRNA_harmony_group2_filtered@meta.data %>%
  dplyr::rename(sample = group2) %>%
  mutate(sample = recode(sample,
                        "Without" = "Without treatment",
                        "Treat" = "With treatment")) %>%
  mutate(sample = factor(sample, levels = c("Without treatment", "With treatment"))) %>%
  group_by(sample, celltype) %>%
  summarize(count = n()) %>%
  spread(celltype, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('sample', 'total_cell_count', everything())) %>%
  arrange(factor(sample, levels = c("Without treatment", "With treatment")))

p3 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable),
           position = 'fill',
           stat = 'identity',
           width = 0.9) +
  scale_fill_manual(name = 'Seurat Clusters', values = custom_colors) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

ggsave(filename = "./data/FigS9/02.cell_count_group2_ALL.pdf", p3, width = 4, height = 8)
ggsave(filename = "./data/FigS9/02.cell_count_group2_ALL.png", p3, width = 4, height = 8, dpi = 600)

# ------------------------------------------------------------------------------
# Proportion by Group (EGPA, Remission, Relapsing)
# ------------------------------------------------------------------------------
table_samples_by_clusters <- scRNA_harmony_group2_filtered@meta.data %>%
  dplyr::rename(sample = group) %>%
  mutate(sample = recode(sample,
                        "E" = "EGPA",
                        "Relapsing" = "Cs-relapsing",
                        "Remission" = "Cs-remission")) %>%
  mutate(sample = factor(sample, levels = c("EGPA", "Cs-remission", "Cs-relapsing"))) %>%
  group_by(sample, celltype) %>%
  summarize(count = n()) %>%
  spread(celltype, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('sample', 'total_cell_count', everything())) %>%
  arrange(factor(sample, levels = levels(scRNA_harmony_group2_filtered@meta.data$sample)))

p6 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable),
           position = 'fill',
           stat = 'identity',
           width = 0.9) +
  scale_fill_manual(name = 'Seurat Clusters', values = custom_colors) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

ggsave(filename = "./data/FigS9/02.cell_count_group_ALL.pdf", p6, width = 5, height = 8)
ggsave(filename = "./data/FigS9/02.cell_count_group_ALL.png", p6, width = 5, height = 8, dpi = 600)

# ==============================================================================
# FIGURE S9G: Relative Enrichment Heatmap (Ro/e)
# ==============================================================================

load("./data/FigS9/scRNA_harmony_group2_filtered.rdata")

# Subset by treatment status
scRNA_harmony_group2_With <- subset(scRNA_harmony_group2_filtered,
                                    subset = group2 == c('With treatment'))
scRNA_harmony_group2_Without <- subset(scRNA_harmony_group2_filtered,
                                       subset = group2 == c('Without treatment'))

# Prepare metadata
treatment <- scRNA_harmony_group2_With@meta.data
treatment <- treatment[, c("group", "celltype")]

Withou <- scRNA_harmony_group2_Without@meta.data
Withou <- Withou[, c("group", "celltype")]

metaTotal <- rbind(Withou, treatment)

# Calculate relative abundance
abundance <- table(metaTotal$celltype, metaTotal$group)

# Calculate total cells and cell type proportions
total_cells <- sum(abundance)
cell_type_proportions <- colSums(abundance) / total_cells

# Calculate expected counts
expected_counts <- outer(cell_type_proportions, rowSums(abundance), "*")

# Calculate observed/expected ratio (Ro/e)
enrichment_ratio <- abundance / t(expected_counts)

RoeResult <- enrichment_ratio[c("Macrophage1", "Macrophage2", "IGF1 Macrophage",
                                "Mitotic macrophage", "mDC1", "mDC2", "XCXR1 DC",
                                "GPR34 DC", "EREG DC", "TREM2 DC", "pDC", "DC (in PBMC)"), ]
RoeResult <- RoeResult[, c("EGPA", "Cs-remission", "Cs-relapsing")]

# Create heatmap
p1 <- pheatmap(RoeResult,
               cluster_rows = F,
               cluster_cols = F,
               color = colorRampPalette(c(brewer.pal(9, "Reds")))(100),
               display_numbers = T,
               angle_col = 45)

pdf(file = "./data/FigS9/07.RoeHeatmap.pdf", width = 3, height = 5.5)
pheatmap(RoeResult,
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c(brewer.pal(9, "Reds")))(100),
         display_numbers = T,
         angle_col = 45)
dev.off()

png(file = "./data/FigS9/07.RoeHeatmap.png", width = 3, height = 5.5)
pheatmap(RoeResult,
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c(brewer.pal(9, "Reds")))(100),
         display_numbers = T,
         angle_col = 45)
dev.off()

# ==============================================================================
# FIGURE S9F: Marker Gene Heatmap
# ==============================================================================

load("./data/FigS9/scRNA_harmony_group2_filtered.rdata")

scRNA_harmony_group2_filtered$celltype <- factor(
  scRNA_harmony_group2_filtered$celltype,
  levels = c("Macrophage1", "Macrophage2", "IGF1 Macrophage", "Mitotic macrophage",
            "mDC1", "mDC2", "XCXR1 DC", "GPR34 DC", "EREG DC", "TREM2 DC",
            "pDC", "DC (in PBMC)")
)

# Find all markers
FindAllMarkers(scRNA_harmony_group2_filtered) -> cell.allmarkers
a <- cell.allmarkers[which(cell.allmarkers[, 2] >= 0.25 & cell.allmarkers[, 5] <= 0.05), ]
a <- as.data.frame(a)

write.table(a, "./data/FigS9/07.findMarker.csv", sep = ',', row.names = F, quote = F)

# Select top 10 markers per cluster
top10 <- a %>%
  group_by(cluster) %>%
  dplyr::slice(1:10)

write.table(top10, "./data/FigS9/07.findMarker_top10.csv", sep = ',', row.names = F, quote = F)

# Create heatmap
p <- DoHeatmap(scRNA_harmony_group2_filtered,
               features = top10$gene,
               disp.min = -2,
               disp.max = 2)

bks <- seq(-2.1, 2.1, by = 0.1)
heatmapColors <- colorRampPalette(rev(c(brewer.pal(7, "RdBu"))))(length(bks) - 1)

p1 <- p + scale_fill_gradientn(colors = heatmapColors)

pdf(file = "./data/FigS9/07.Marker_Heatmap.pdf", width = 8, height = 6)
p + scale_fill_gradientn(colors = heatmapColors) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )
dev.off()

# ------------------------------------------------------------------------------
# Heatmap with Custom Gene List
# ------------------------------------------------------------------------------
gene <- c("FN1", "MRC1", "MARCO", "CD9", "OLR1", "SCGB1A1", "FOLR3", "SLPI", "FCGR3B",
          "TOP2A", "MKI67", "CCL22", "CCR7", "LAMP3", "GPR183", "CCL17", "CD1C", "CD1E",
          "FCER1A", "CD207", "CLEC10A", "XCR1", "WDFY4", "CLEC9A", "CXCL9", "CXCL10",
          "TMEM176B", "TMEM176A", "SPP1", "CHI3L1", "MMP9", "IL1RN", "MMP12", "CCL2",
          "CXCL8", "CCL3", "VEGFA", "TREM2", "PLD4", "IRF7", "IRF8", "FOS", "JUNB", "CXCR4")

# Read top genes from Excel
top <- read_excel("./data/FigS9/07.top.xlsx")
gene1 <- c(top$gene, gene)

p1 <- DoHeatmap(scRNA_harmony_group2_filtered,
                features = top$gene,
                disp.min = -2,
                disp.max = 2)

pdf(file = "./data/FigS9/07.Marker_Heatmap_gene.pdf", width = 8, height = 6)
p1 + scale_fill_gradientn(colors = heatmapColors) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )
dev.off()

# ==============================================================================
# End of Script
# ==============================================================================
