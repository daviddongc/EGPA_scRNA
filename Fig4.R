# ==============================================================================
# EGPA scRNA-seq Analysis: Figure 4
# ==============================================================================
# Description: This script generates Figure 4
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic 1
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
macroTreat <- readRDS("./data/Fig4/macroTreat.RDS")

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
macro <- readRDS("./data/Fig4/macro.RDS")
macro$group2 <- "Without"

# Merge baseline and treatment datasets
scRNA <- merge(macro, macroTreat)
system.time(save(scRNA, file = "./data/Fig4/scRNA.orig.rdata"))

# ------------------------------------------------------------------------------
# Process Combined Dataset
# ------------------------------------------------------------------------------
rm(list = ls())
load("./data/Fig4/scRNA.orig.rdata")

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


# ------------------------------------------------------------------------------
# Harmony Batch Correction
# -----------------------------------------------------------------------------

# Harmony correction by group2 (treatment status)
scRNA_harmony_group2 <- RunHarmony(scRNA_ERR,
                                   reduction = "pca",
                                   group.by.vars = "group2",
                                   reduction.save = "harmony",
                                   num.cores = 55)
FindNeighbors(scRNA_harmony_group2, dims = c(1:19), reduction = "harmony") -> scRNA_harmony_group2
FindClusters(scRNA_harmony_group2, resolution = 1, reduction = "harmony") -> scRNA_harmony_group2
RunUMAP(scRNA_harmony_group2, dims = c(1:19), reduction = "harmony") -> scRNA_harmony_group2


# ------------------------------------------------------------------------------
# Differential Expression Analysis
# ------------------------------------------------------------------------------

# Find markers using presto
table(scRNA_harmony_group2$SCT_snn_res.1)
system.time(dge <- wilcoxauc(scRNA_harmony_group2, group_by = 'SCT_snn_res.1'))
dge_0.5 <- dge[which(dge$logFC > 0.05 | dge$logFC < -0.5), ]
write.table(dge_0.5, file = './data/Fig4/AllMarkers_presto.csv', sep = ',', row.names = FALSE, quote = FALSE)

# Find markers using Seurat
Idents(scRNA_harmony_group2) <- scRNA_harmony_group2$SCT_snn_res.1
system.time(dge1 <- FindAllMarkers(scRNA_harmony_group2, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = 0.25, test.use = 'wilcox',
                                   return.thresh = 0.01))
write.table(dge1, file = './data/Fig4/AllMarkers_findmarker.csv', sep = ',', row.names = FALSE, quote = FALSE)

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
system.time(save(scRNA_harmony_group2_filtered, file = "./data/Fig4/scRNA_harmony_group2_filtered.rdata"))
# ==============================================================================
# FIGURE 4A: UMAP Plot Split by Treatment Status
# ==============================================================================
custom_colors <- c("Macrophage1"="#A6CEE3","Macrophage2"="#FB9A99","IGF1 Macrophage"="#1F78B4",
                   "Mitotic macrophage"="#E31A1C","mDC1"="#CAB2D6","mDC2"="darkred","XCXR1 DC"="#33A02C",
                   "GPR34 DC"="#FDBF6F","EREG DC"="darkgrey","TREM2 DC"="#7570B3","pDC"="#D95F02","DC (in PBMC)"="#6A3D9A","#B15928","#1B9E77","#FF7F00","#FFFF99","black","navy","firebrick")

anno_plot2 <- DimPlot2(scRNA_harmony_group2_filtered,
                       reduction = "umap",
                       split.by = "group2",
                       label = F,
                       label.size = 2.8,
                       cols = custom_colors,
                       theme = theme_umap_arrows()) +
  labs(x = "UMAP1", y = "UMAP2", title = "UMAP") +
  theme(
    axis.title.x = element_text(size = 20, family = "Times", face = "bold"),
    axis.text.x = element_text(angle = 0, size = 20, hjust = 0.5, family = "Times", face = "bold"),
    axis.title.y = element_text(size = 20, family = "Times", face = "bold"),
    axis.text.y = element_text(size = 16, family = "Times", face = "bold"),
    legend.position = "none"
  )

ggsave(filename = './data/Fig4/UMAP_group2.pdf', anno_plot2, w = 14, h = 8)
ggsave(filename = './data/Fig4/UMAP_group2.png', anno_plot2, w = 14, h = 8)

# ==============================================================================
# FIGURE 4B: Volcano Plot for IGF1 Macrophage
# ==============================================================================

# ------------------------------------------------------------------------------
# Find Markers in Treatment Group
# ------------------------------------------------------------------------------
scRNA_harmony_group2_Treat <- subset(scRNA_harmony_group2_filtered, subset = group2 == c('With treatment'))
Idents(scRNA_harmony_group2_Treat) <- scRNA_harmony_group2_Treat$celltype
system.time(Treat <- FindAllMarkers(scRNA_harmony_group2_Treat,
                                    logfc.threshold = 0.25,
                                    min.pct = 0.25,
                                    test.use = 'wilcox',
                                    return.thresh = 0.01))
write.table(Treat, file = './data/Fig4/Treat_findallmarkers.csv', sep = ',', row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# Prepare Data for Volcano Plot
# ------------------------------------------------------------------------------
DEGs <- Treat %>%
  dplyr::filter(., cluster == "IGF1 Macrophage") %>%
  mutate(., change = case_when(avg_log2FC >= .5 ~ "up",
                               avg_log2FC <= -0.5 ~ "down",
                               TRUE ~ "normal"),
         symbol = gene) %>%
  dplyr::rename(group = cluster,
                log2FoldChange = avg_log2FC)

# ------------------------------------------------------------------------------
# Create Volcano Plot
# ------------------------------------------------------------------------------
p <- ggplot(data = DEGs) +
  geom_jitter(data = DEGs %>%
                dplyr::filter(change == "normal"),
              aes(x = group, y = log2FoldChange,
                  color = change, size = abs(log2FoldChange),
                  alpha = abs(log2FoldChange)),
              width = 0.4) +
  geom_jitter(data = DEGs %>%
                dplyr::filter(change != "normal"),
              aes(x = group, y = log2FoldChange,
                  color = change, size = abs(log2FoldChange),
                  alpha = abs(log2FoldChange)),
              width = 0.4) +
  geom_text_repel(data = DEGs %>%
                    dplyr::filter(change == "up") %>%
                    dplyr::group_by(group) %>%
                    dplyr::arrange(desc(abs(log2FoldChange))) %>%
                    dplyr::slice_head(n = 10) %>%
                    dplyr::ungroup() %>%
                    dplyr::distinct(),
                  aes(x = group, y = log2FoldChange, label = symbol)) +
  geom_tile(aes(x = group, y = 0, fill = group), height = 0.5) +
  geom_text(data = DEGs %>%
              dplyr::select(group) %>%
              dplyr::distinct(group, .keep_all = T),
            aes(x = group, y = 0, label = group),
            size = 4) +
  scale_y_continuous(limits = c(-2.5, 3.1)) +
  scale_size(range = c(.5, 8)) +
  scale_alpha(range = c(0.1, 1)) +
  scale_color_manual(values = c("up" = "#f46d43",
                                "normal" = "#bdbdbd",
                                "down" = "#3288bd")) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Cell Type") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "#000000", size = 12),
    axis.title = element_text(color = "#000000", size = 15),
    panel.grid = element_blank(),
    legend.background = element_roundrect(color = "#969696")
  )

ggsave(filename = "./data/Fig4/Treat_Volcano_IGF1.pdf", p, width = 6, height = 4)
ggsave(filename = "./data/Fig4/Treat_Volcano_IGF1.png", p, width = 6, height = 4, dpi = 600)

# ==============================================================================
# FIGURE 4C: GO Enrichment Analysis
# ==============================================================================
#FindMarker & GO
Idents(scRNA_harmony_group2_Treat) <- scRNA_harmony_group2_Treat$celltype
system.time(GOgene <- FindMarkers(scRNA_harmony_group2_Treat,ident.1 = "IGF1 Macrophage",only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25,test.use = 'wilcox',return.thresh = 0.01))
write.table(GOgene,file = './data/Fig4/IGF1_findmarker.csv',sep = ',',row.names = FALSE,quote = FALSE)

# Query available databases
dbs <- enrichR::listEnrichrDbs()
print(dbs)

# Prepare gene list (top 20 upregulated genes)
gene <- unique(rownames(GOgene %>%
                          dplyr::arrange(desc(avg_log2FC)) %>%
                          dplyr::slice_head(n = 20)))

gene_list <- gene

# Select database
selected_dbs <- c("Reactome_Pathways_2024")

# Run enrichment analysis
results <- enrichr(gene_list, databases = selected_dbs)
print(results)

# Extract results
Reactome <- results$Reactome_Pathways_2024
write.table(Reactome, file = './data/Fig4/Reactome_Pathways_2024.csv', sep = ',', row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# Create Enrichment Plot
# ------------------------------------------------------------------------------
GO <- Reactome %>%
  arrange(P.value) %>%
  slice_head(n = 5) %>%
  mutate(x = seq(1:5)) %>%
  mutate(Pvalue = -log10(P.value))

p <- ggplot(GO, aes(x = Pvalue, y = reorder(Term, Pvalue), fill = Pvalue)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lightgrey", high = "#1F78B4") +
  geom_text(aes(x = 0.01, label = Term), hjust = 0, fontface = "bold",
            family = "Times", color = "black", size = 5.5) +
  geom_text(aes(x = 0.01, y = Term, label = Genes),
            size = 4, fontface = "italic", hjust = 0, vjust = 5,
            family = "Times", color = "black", lineheight = 0.9) +
  labs(x = "-Log10 P Value",
       y = "KEGG Pathways",
       title = "TOP20 genes overexpressed by IGF1+ macrophages",
       fill = "-log10(p-value)") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 14),
    axis.title.x = element_text(face = "bold", color = "black", family = "Times", size = 18),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12),
    legend.title = element_text(face = "bold", color = "black", family = "Times", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 20),
    axis.ticks.y = element_blank()
  ) +
  scale_y_discrete(expand = c(0.05, 0)) +
  guides(fill = "none")

ggsave(filename = "./data/Fig4/Reactome_Top20.pdf", p, width = 8, height = 6)
ggsave(filename = "./data/Fig4/Reactome_Top20.png", p, width = 8, height = 6, dpi = 600)

# ==============================================================================
# FIGURE 4D: Feature Plots for IGF1, HP, RETN
# ==============================================================================

# Subset to treatment group
scRNA_harmony_group2_Treat <- subset(scRNA_harmony_group2_filtered,
                                     subset = group2 == c('With treatment'))

# IGF1 expression
p <- FeaturePlot(scRNA_harmony_group2_Treat,
                 features = c("IGF1"),
                 min.cutoff = 1,
                 max.cutoff = 7) &
  scale_color_gradientn(
    colours = c('#E2F3F6', '#EFF8DF', '#FEFFC1', '#FFEFA7', '#FEE192', '#FEC97A',
                '#FDB062', '#FA9053', '#F46D43', '#D83328', '#BF1D27', '#AE0C26'),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = NA)
  )

ggsave(filename = './data/Fig4/IGF1_scatter.pdf', p, w = 6.5, h = 6)
ggsave(filename = './data/Fig4/IGF1_scatter.png', p, w = 6.5, h = 6, dpi = 600)

# HP expression
p <- FeaturePlot(scRNA_harmony_group2_Treat,
                 features = c("HP"),
                 min.cutoff = 1,
                 max.cutoff = 7) &
  scale_color_gradientn(
    colours = c('#E2F3F6', '#EFF8DF', '#FEFFC1', '#FFEFA7', '#FEE192', '#FEC97A',
                '#FDB062', '#FA9053', '#F46D43', '#D83328', '#BF1D27', '#AE0C26'),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = NA)
  )

ggsave(filename = './data/Fig4/HP_scatter.pdf', p, w = 6.5, h = 6)
ggsave(filename = './data/Fig4/HP_scatter.png', p, w = 6.5, h = 6, dpi = 600)

# RETN expression
p <- FeaturePlot(scRNA_harmony_group2_Treat,
                 features = c("RETN"),
                 min.cutoff = 1,
                 max.cutoff = 7) &
  scale_color_gradientn(
    colours = c('#E2F3F6', '#EFF8DF', '#FEFFC1', '#FFEFA7', '#FEE192', '#FEC97A',
                '#FDB062', '#FA9053', '#F46D43', '#D83328', '#BF1D27', '#AE0C26'),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = NA)
  )

ggsave(filename = './data/Fig4/RETN_scatter.pdf', p, w = 6.5, h = 6)
ggsave(filename = './data/Fig4/RETN_scatter.png', p, w = 6.5, h = 6, dpi = 600)

# ==============================================================================
# FIGURE 4E: Violin Plots with MAGIC Imputation
# ==============================================================================

scRNA_harmony_group2_filtered <- readRDS("./data/Fig4/scRNA_harmony_group2_filtered.RDS")

# Extract expression matrix
expression_matrix <- as.matrix(scRNA_harmony_group2_filtered[["SCT"]]@data)

# ------------------------------------------------------------------------------
# Configure Python Environment for MAGIC
# ------------------------------------------------------------------------------
# NOTE: Configure your own conda environment and Python path before running MAGIC
# Set conda path
# conda_path <- "/path/to/your/conda"
# env_path <- "/path/to/your/magic/env"
#
# # Activate conda environment
# use_condaenv(env_path, conda = conda_path, required = TRUE)
# use_python("/path/to/your/python", required = TRUE)
# py_config()

# ------------------------------------------------------------------------------
# Run MAGIC Imputation
# ------------------------------------------------------------------------------
MAGIC_data <- magic(expression_matrix)

# Add MAGIC results to Seurat object
scRNA_harmony_group2_filtered[["MAGIC_STC"]] <- CreateAssayObject(MAGIC_data$result)

# ------------------------------------------------------------------------------
# Violin Plot for pDC
# ------------------------------------------------------------------------------
spe <- subset(scRNA_harmony_group2_filtered, idents = "pDC")
Idents(spe) <- "group"
spe@active.assay <- "MAGIC_STC"

plots <- VlnPlot(spe,
                 features = c("IRF7", "IRF8", "MX1", "PLD4"),
                 pt.size = 0,
                 combine = FALSE,
                 adjust = 1)
CombinePlots(plots = plots, ncol = 1)

p3 <- VlnPlot(spe,
              stack = TRUE,
              features = c("IRF7", "IRF8", "MX1", "PLD4"),
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

# ------------------------------------------------------------------------------
# Violin Plot for Macrophage
# ------------------------------------------------------------------------------
Idents(scRNA_harmony_group2_filtered) <- "celltype"
Macrophage <- subset(scRNA_harmony_group2_filtered, idents = "Macrophage")
Idents(Macrophage) <- "group"
Macrophage@active.assay <- "MAGIC_STC"

plots <- VlnPlot(Macrophage,
                 features = c("HP", "IGF1", "RETN"),
                 pt.size = 0,
                 combine = FALSE,
                 adjust = 1,
                 cols = c("#BC1A2A", "#F09C00", "#EF90B1")) +
  theme(legend.position = "none")
CombinePlots(plots = plots, ncol = 1)


# clear environment
rm(list = ls())

# load necessary libraries
library(ggplot2)
library(Seurat)
library(sctransform)
library(RColorBrewer)
library(patchwork)
library(pheatmap)
library(harmony)
library(dplyr)
library(future)
library(furrr)
library(presto)
library(pheatmap)
library(reticulate)

# ==============================================================================
# FIGURE 4G: UMAP Split by Treatment
# ==============================================================================
#### Figure4g
scRNA_harmony_group_filtered <- readRDS("./data/Fig4/F4_g.rds")
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

p2<-DimPlot(scRNA_harmony_group_filtered,
        label =F,
        label.size = 2.5,
        group.by = "cell_type",  # Group by cell_type
        split.by = "treatment",  # Split by treatment condition
        pt.size = 0.3) +  # Point size for cells
  scale_color_manual(values = custom_palette) +  # Apply custom color palette
  labs(color = "Cell Type") +  # Set legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
p2
ggsave(filename = './data/Fig4/figure4_g.png', p2, w = 5.5, h = 3)
ggsave(filename = './data/Fig4/figure4_g.pdf', p2, w = 5.5, h = 3)

# ==============================================================================
# FIGURE 4I: Dotplot of Hormones and Receptors
# ==============================================================================
#### Figure4i
epi <- readRDS("./data/Fig4/F4_i.rds")
epi$celltype<-epi$cell_type

hormones<-c("SST","CALCA","AGT","IGF1","NPPC","EDN1")
horReceptors<-c("OXTR","EDNRA","IGF1R","PTGER4","VIPR1","PTGER3","ADRB2")

## factor----levels
epi@active.ident<-factor(epi$celltype,levels=rev(c("Basal",
                                                   "Club cells",
                                                   "Cycling",
                                                   "Deuterosome-stage ciliated",
                                                   "Differentiating basal",
                                                   'Goblet-1',
                                                   "Goblet-2",
                                                   "Ionocyte",
                                                   "Mature ciliated",
                                                   "Neuroendocrine",
                                                   "Proximal ciliated")))
p3<-DotPlot(epi,features=horReceptors,
            cols=c("lightgrey","#252525"),
            col.min=-1,col.max=1,
            dot.min=0,dot.scale=6,
            scale.min=0,scale.max=50)+
  theme(
    axis.text.x = element_text(vjust = 0.5, angle = -90), ##   hjust = 0.5, 
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.line = element_line(color = "black", size = 0.5), 
    plot.margin = margin(10, 10, 10, 10)  
  ) ;p3

p4<-DotPlot(epi,features=hormones,
            cols=c("lightgrey","#252525"),
            col.min=-1,col.max=1, 
            dot.min=0,dot.scale=6,
            scale.min=0,scale.max=50)+
  theme(
    axis.text.x = element_text(vjust = 0.5, angle = -90), ##   hjust = 0.5, 
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.line = element_line(color = "black", size = 0.5), 
    plot.margin = margin(10, 10, 10, 10)
  ) ;p4

ggsave(filename = './data/Fig4/F4_i_horReceptors.pdf',p3,w=6,h=4)
ggsave(filename = './data/Fig4/F4_i_horReceptors.png',p3,w=6,h=4,dpi = 300)
ggsave(filename = './data/Fig4/F4_i_hormones.pdf',p4,w=6,h=4)
ggsave(filename = './data/Fig4/F4_i_hormones.png',p4,w=6,h=4,dpi = 300)

# ==============================================================================
# FIGURE 4J: Violin Plots for Goblet and Ionocyte
# ==============================================================================
#### Figure4j
MAGIC_data <- readRDS("./data/Fig4/F4_j.rds")
#### 01 Goblet_1 ####
Idents(MAGIC_data)<-"celltype"
Goblet_1<-subset(MAGIC_data,idents=c("Goblet-1"))  ## Goblet_1   Goblet_2    lonocyte
Goblet_1@active.assay<-"MAGIC_SCT"  ## MAGIC_SCT
table(Goblet_1$celltype)
Idents(Goblet_1)<-'group'
table(Goblet_1$group)
Goblet_1$group <- factor(Goblet_1$group, levels = c("E",'Remission',"Relapsing"))
Idents(Goblet_1)<- "group"

p111<-VlnPlot(Goblet_1,
              stack=TRUE,
              features=c('MUC2','MUC5AC','IGFBP5','IGFBP7'),
              ident = c("E",'Remission',"Relapsing"),
              flip=T,
              same.y.lims=FALSE,
              fill.by='ident',
              cols=c('#b10025','#f48515','#ef719f'))+
  scale_y_continuous(limits = c(0, 2)) +
  ggtitle("Goblet_1") + 
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank()
  );p111


#### 02 Goblet_2 ####
Idents(MAGIC_data)<-"celltype"
Goblet_2<-subset(MAGIC_data,idents=c("Goblet-2"))  ## Goblet_1   Goblet_2    lonocyte
Goblet_2@active.assay<-"MAGIC_SCT"  ## RNA   MAGIC_SCT
table(Goblet_2$celltype)
Idents(Goblet_2)<-'group'
table(Goblet_2$group)
Goblet_2$group <- factor(Goblet_2$group, levels = c("E",'Remission',"Relapsing"))
Idents(Goblet_2)<- "group"

p22<-VlnPlot(Goblet_2,
             stack=TRUE,
             features=c('IDO1','CXCL9','CXCL10','IL19'),
             ident = c("E",'Remission',"Relapsing"),
             flip=T,
             same.y.lims=FALSE,
             fill.by='ident',
             cols=c('#b10025','#f48515','#ef719f'))+

  scale_y_continuous(limits = c(0, 0.1)) +
  ggtitle("Goblet_2") +  # 添加标题
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank()
  );p22

####  lonocyte  ####
Idents(MAGIC_data)<-"celltype"
lonocyte<-subset(MAGIC_data,idents=c("Ionocyte"))  ## Goblet_1   Goblet_2    lonocyte
lonocyte@active.assay<-"MAGIC_SCT"   ##  MAGIC_SCT
table(lonocyte$celltype)
Idents(lonocyte)<-'group'
table(lonocyte$group)
lonocyte$group <- factor(lonocyte$group, levels = c("E",'Remission',"Relapsing"))
Idents(lonocyte)<- "group"
p3<-VlnPlot(lonocyte,
            stack=TRUE,
            features=c('CFTR','IGF1'),
            ident = c("E",'Remission',"Relapsing"),
            flip=T,
            same.y.lims=FALSE,
            fill.by='ident',
            cols=c('#b10025','#f48515','#ef719f'))+

  scale_y_continuous(limits = c(0,2.2)) +
  ggtitle("lonocyte") +  # 添加标题
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank());p3

ggsave(filename = './data/Fig4/Goblet_1.pdf',p111,w=3,h=4)
ggsave(filename = './data/Fig4/Goblet_1.png',p111,w=3,h=4,dpi = 600)
ggsave(filename = './data/Fig4/Goblet_2.pdf',p22,w=3,h=4)
ggsave(filename = './data/Fig4/Goblet_2.png',p22,w=3,h=4,dpi = 600)
ggsave(filename = './data/Fig4/Ionocyte.pdf',p3,w=3,h=3)
ggsave(filename = './data/Fig4/Ionocyte.png',p3,w=3,h=3,dpi = 600)

# ==============================================================================
# FIGURE 4H: Heatmap of Cell Type Enrichment
# ==============================================================================
#### Figure4h
scRNA_harmony_group2_filtered <- readRDS("./data/Fig4/F4_i.rds")
scRNA_harmony_group2_filtered$celltype<-scRNA_harmony_group2_filtered$cell_type

scRNA_harmony_group2_filtered$Type <- ifelse(
  is.na(scRNA_harmony_group2_filtered$type), 
  scRNA_harmony_group2_filtered$source, 
  scRNA_harmony_group2_filtered$type
)

##  extract Biopsy
Idents(scRNA_harmony_group2_filtered)<-'Type'
scRNA_harmony_group2_filtered<-subset(scRNA_harmony_group2_filtered,idents= c('Biopsy'))

## rename the group
scRNA_harmony_group2_filtered$group <- ifelse(scRNA_harmony_group2_filtered$group == 'E', 'EGPA',
                                              ifelse(scRNA_harmony_group2_filtered$group == 'Cs-relapsing', 'Cs-remission', 
                                                     scRNA_harmony_group2_filtered$group))
scRNA_harmony_group2_filtered$group <- factor(scRNA_harmony_group2_filtered$group, levels = c("EGPA","Relapsing",'Remission'))


scRNA_harmony_group2_With <- subset(scRNA_harmony_group2_filtered, subset = treatment == c('With_treatment'))
treatment<-scRNA_harmony_group2_With@meta.data
treatment<-treatment[,c("group","celltype")];head(treatment,2)

scRNA_harmony_group2_Without <- subset(scRNA_harmony_group2_filtered, subset = treatment == c('without_treatment'))
Withou<-scRNA_harmony_group2_Without@meta.data
Withou<-Withou[,c("group","celltype")];head(Withou,2)


metaTotal<-rbind(Withou,treatment);head(metaTotal,2)

abundance <- table(metaTotal$celltype, metaTotal$group);head(abundance,2)
total_cells <- sum(abundance);total_cells
cell_type_proportions <- colSums(abundance) / total_cells;cell_type_proportions
expected_counts <- outer(cell_type_proportions, rowSums(abundance), "*") ;expected_counts
enrichment_ratio <- abundance / t(expected_counts) ;enrichment_ratio
RoeResult<-enrichment_ratio[c("Basal","Club cells","Cycling",
                              "Deuterosome-stage ciliated",'Differentiating basal',
                              "Goblet-1","Goblet-2","Ionocyte",
                              "Mature ciliated","Neuroendocrine",'Proximal ciliated'),]
RoeResult<-RoeResult[,c("EGPA","Remission","Relapsing")]
colnames(RoeResult)
colnames(RoeResult)<-c("EGPA","Cs-remission","Cs-relapsing")
RoeResult <- RoeResult[, c("Cs-remission","Cs-relapsing",  "EGPA")]

RoeResult_1<-RoeResult[,c(2,1,3)]
p1<-pheatmap(RoeResult_1,
             cluster_rows=F,cluster_cols=F,
             cellwidth = 25, cellheight = 25,
             color=colorRampPalette(c(brewer.pal(9,"Reds")))(100),
             display_numbers=T);p1
ggsave(filename = './data/Fig4/Figure4h.pdf',p1,w=4,h=6)
ggsave(filename = './data/Fig4/Figure4h.png',p1,w=4,h=6,dpi = 300)

# ==============================================================================
# End of Script
# ==============================================================================
