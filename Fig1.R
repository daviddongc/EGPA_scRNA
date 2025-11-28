# ==============================================================================
# EGPA scRNA-seq Analysis: Figure 1 
# ==============================================================================
# Description: This script generates Figure 1C, 1D 
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
# library(xlsx)
library(dplyr)
library(tidyr)
library(rstatix)
library(FSA)  # For Dunn's test
library(ggplot2)
library(ggsignif)
library(ggpubr)

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
combined <- readRDS("./data/Fig1/combined.RDS")

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
# FIGURE 1C: UMAP Plot of All Cell Types
# ==============================================================================
p <- DimPlot(combined, reduction = "umap", cols = mypalette, pt.size = 1, label = FALSE)

# ==============================================================================
# FIGURE 1D: Cell Type Composition Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# Prepare Proportion Data
# ------------------------------------------------------------------------------
pdata <- data.frame(table(combined$subject_id, combined$source, combined@active.ident))

# Add group information
pdata$Group <- rep(0, length(pdata[, 1]))
for (i in 1:length(pdata[, 1])) {
  pdata[i, "Group"] <- as.character(unique(combined@meta.data[which(combined$subject_id == pdata[i, 1]), "group"]))
}

# Calculate total cells per subject and source
pdata$Total <- rep(0, length(pdata[, 1]))
for (i in 1:length(pdata[, 1])) {
  pdata[i, "Total"] <- length(combined@meta.data[which(combined$subject_id == pdata[i, 1] &
                                                         combined$source == pdata[i, 2]), 1])
}
pdata$Proportion <- pdata$Freq / pdata$Total

# Set subject order
pdata$Var1 <- factor(pdata$Var1,
                     levels = c("P1-ZJ", "P2-XFT", "P3-ZYL", "P4-CBX", "P5-YQJ", "P6-RL",
                               "P1-LL", "P2-LZQ", "P3-HXQ", "P4-FKY", "P5-CZH", "P6-ZSQ",
                               "P7-LBZ", "P8-YL", "P9-ZXY", "P10-GQ", "P11-LLX", "P12-HYF", "A-1",
                               "P1-ZQJ", "P2-HYC", "P3-LQX", "P4-ZYX", "P5-WYX", "P6-LHQ",
                               "P7-JK", "P8-WL", "E-1"))

# Rename "Mono DCs" to "DC" in factor labels
pdata$Var3 <- factor(pdata$Var3,
                     levels = c("Epithelial", "Neuroendocrine", "Fibroblast",
                               "Smooth muscle", "Endothelial", "T & NK cells",
                               "B cells", "Myeloblast", "Macrophage", "Mono DCs",
                               "pDC", "Neutrophils-1", "Neutrophils-2",
                               "Eosinophils", "Mast cells", "Megakaryocyte",
                               "Erythroid cells"),
                     labels = c("Epithelial", "Neuroendocrine", "Fibroblast",
                               "Smooth muscle", "Endothelial", "T & NK cells",
                               "B cells", "Myeloblast", "Macrophage", "DC",
                               "pDC", "Neutrophils-1", "Neutrophils-2",
                               "Eosinophils", "Mast cells", "Megakaryocyte",
                               "Erythroid cells"))

# Save proportion data
saveRDS(pdata, "./data/Fig1/pdata.RDS")

# ------------------------------------------------------------------------------
# BALF Analysis - Stacked Bar Plot
# ------------------------------------------------------------------------------
tissue <- "BALF"
data <- pdata[which(pdata$Var2 == tissue), ]
data <- data[which(!is.na(data$Proportion)), ]

# Calculate normalized proportions per group
data_summary <- data %>%
  group_by(Group, Var3) %>%
  summarise(Proportion = sum(Proportion, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Proportion = Proportion / sum(Proportion)) %>%
  ungroup()

# Recode group labels
data_summary$Group <- factor(data_summary$Group,
                             levels = c("E", "A", "C"),
                             labels = c("EGPA", "Asthma", "Control"))

# Create stacked bar plot
ggplot(data_summary, aes(x = Group, y = Proportion, fill = Var3)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_family = "Arial") +
  scale_fill_manual(values = mypalette) +
  labs(x = "", y = "Proportion", fill = "") +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", family = "Arial", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14, face = "bold", family = "Arial"),
    axis.title = element_text(size = 16, face = "bold", family = "Arial"),
    legend.text = element_text(size = 14, face = "bold", family = "Arial"),
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),
    panel.border = element_rect(color = "black", size = 1.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.ticks.length = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(ncol = 6, keywidth = 0.8, keyheight = 0.8)) +
  coord_flip()

# ------------------------------------------------------------------------------
# BALF Analysis - Box Plots with Statistics
# ------------------------------------------------------------------------------

# Recode group labels
# Prepare BALF data
tissue <- "BALF"
data <- pdata[which(pdata$Var2 == tissue), ]
data <- data[which(!is.na(data$Proportion)), ]

# Recode group labels
data <- data %>%
  mutate(Group = recode(Group,
                        "C" = "Control",
                        "A" = "Asthma",
                        "E" = "EGPA")) %>%
  dplyr::rename(CellType = Var3)

# Set factor levels
data$Group <- factor(data$Group, levels = c("Control", "Asthma", "EGPA"))

# ------------------------------------------------------------------------------
# Kruskal-Wallis Test with Effect Size (η²H)
# ------------------------------------------------------------------------------

# Function to calculate Kruskal-Wallis with effect size
perform_kw_test <- function(data) {
  # Get valid cell types (present in all 3 groups)
  valid_celltypes <- data %>%
    filter(Proportion > 0) %>%
    group_by(CellType) %>%
    summarise(n_groups = n_distinct(Group)) %>%
    filter(n_groups >= 3) %>%
    pull(CellType)
  
  # Filter data for valid cell types
  filtered_data <- data %>%
    filter(CellType %in% valid_celltypes)
  
  # Perform Kruskal-Wallis test for each cell type
  kw_results <- filtered_data %>%
    group_by(CellType) %>%
    kruskal_test(Proportion ~ Group) %>%
    mutate(
      # Calculate eta-squared H (effect size)
      eta_squared_H = statistic / (n - 1)
    )
  
  return(kw_results)
}

kw_results <- perform_kw_test(data)

# Print Kruskal-Wallis results
cat("\n=== Kruskal-Wallis Test Results ===\n")
print(kw_results)

# ------------------------------------------------------------------------------
# Dunn's Post-hoc Test with Bonferroni Adjustment
# ------------------------------------------------------------------------------

# Function to perform Dunn's test
perform_dunn_test <- function(data) {
  # Get valid cell types
  valid_celltypes <- data %>%
    filter(Proportion > 0) %>%
    group_by(CellType) %>%
    summarise(n_groups = n_distinct(Group)) %>%
    filter(n_groups >= 3) %>%
    pull(CellType)
  
  # Filter data for valid cell types
  filtered_data <- data %>%
    filter(CellType %in% valid_celltypes)
  
  # Perform Dunn's test for each cell type
  dunn_results_list <- list()
  
  for (cell_type in valid_celltypes) {
    cell_data <- filtered_data %>%
      filter(CellType == cell_type)
    
    # Perform Dunn's test using FSA package
    dunn_test <- dunnTest(Proportion ~ Group, 
                          data = cell_data, 
                          method = "bh")
    
    # Extract results
    dunn_df <- dunn_test$res %>%
      mutate(CellType = cell_type) %>%
      select(CellType, Comparison, Z, P.unadj, P.adj)
    
    dunn_results_list[[cell_type]] <- dunn_df
  }
  
  # Combine all results
  dunn_results <- bind_rows(dunn_results_list)
  
  # Separate comparison column into group1 and group2
  dunn_results <- dunn_results %>%
    separate(Comparison, into = c("group1", "group2"), sep = " - ", remove = FALSE)
  
  return(dunn_results)
}

dunn_results <- perform_dunn_test(data)

# Print Dunn's test results
cat("\n=== Dunn's Post-hoc Test Results (Bonferroni Adjusted) ===\n")
print(dunn_results)

# ------------------------------------------------------------------------------
# Create Comprehensive Results Table
# ------------------------------------------------------------------------------

# Merge KW and Dunn's results
# First convert kw_results to a regular data frame
kw_df <- as.data.frame(kw_results) %>%
  select(CellType, statistic, df, p, n, eta_squared_H) %>%
  dplyr::rename(KW_chi_squared = statistic,
                KW_df = df,
                KW_p_value = p,
                KW_n = n)

comprehensive_results <- dunn_results %>%
  left_join(kw_df, by = "CellType") %>%
  mutate(
    # Add significance levels
    significance = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01 ~ "**",
      P.adj < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    # Format p-values
    P.adj_formatted = case_when(
      P.adj < 0.001 ~ "< 0.001",
      TRUE ~ sprintf("%.4f", P.adj)
    )
  ) %>%
  arrange(CellType, group1, group2)

# Print comprehensive results
cat("\n=== Comprehensive Statistical Results ===\n")
print(comprehensive_results)

# ------------------------------------------------------------------------------
# Save Results
# ------------------------------------------------------------------------------

# Save Kruskal-Wallis results
write.csv(kw_results, 
          "./BALF_KW_test_results.csv", 
          row.names = FALSE)

# Save Dunn's test results
write.csv(dunn_results, 
          "./BALF_Dunn_test_results.csv", 
          row.names = FALSE)

# Save comprehensive results (for Supplementary Table 5)
write.csv(comprehensive_results, 
          "./BALF_Supplementary_Table5_statistics.csv", 
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Create Publication-Ready Table
# ------------------------------------------------------------------------------

# Format for publication
pub_table <- comprehensive_results %>%
  select(CellType, Comparison, Z, P.adj, P.adj_formatted, significance, 
         KW_chi_squared, KW_p_value, eta_squared_H) %>%
  mutate(
    # Round values for presentation
    Z = round(Z, 3),
    KW_chi_squared = round(KW_chi_squared, 3),
    eta_squared_H = round(eta_squared_H, 3),
    KW_p_value = case_when(
      KW_p_value < 0.001 ~ "< 0.001",
      TRUE ~ sprintf("%.4f", KW_p_value)
    )
  )

# Save publication-ready table
write.csv(pub_table, 
          "./BALF_Publication_Ready_Statistics.csv", 
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Create Summary Report for Figure Legend
# ------------------------------------------------------------------------------

# Extract significant results for figure legend
significant_results <- comprehensive_results %>%
  filter(P.adj < 0.05) %>%
  select(CellType, group1, group2, P.adj_formatted, significance)

# Format for figure legend text
figure_legend_text <- significant_results %>%
  mutate(
    text = sprintf("%s: %s vs %s, P = %s", 
                   CellType, group1, group2, P.adj_formatted)
  ) %>%
  pull(text)

# Save figure legend text
writeLines(figure_legend_text, "./BALF_Figure_Legend_P_values.txt")

cat("\n=== Figure Legend P-values (P < 0.05) ===\n")
for (text in figure_legend_text) {
  cat(paste0("• ", text, "\n"))
}

# ------------------------------------------------------------------------------
# Create Enhanced Boxplot with Dunn's Test Results
# ------------------------------------------------------------------------------

# Filter for selected cell types (in the desired order)
selected_celltypes <- c("Mast cells", "Epithelial", "T & NK cells", 
                        "Neutrophils-2", "Eosinophils", "Macrophage")

filtered_data <- data %>%
  filter(CellType %in% selected_celltypes)

# Function to set GroupCellType order (using custom order)
set_group_celltype_order <- function(data, cell_type_order) {
  groups <- c("Control", "Asthma", "EGPA")
  sorted_levels <- c()
  for (cell_type in cell_type_order) {
    for (group in groups) {
      sorted_levels <- c(sorted_levels, paste(cell_type, group, sep = "."))
    }
  }
  data$GroupCellType <- factor(data$GroupCellType, levels = sorted_levels)
  return(data)
}

# Create GroupCellType column and set order
filtered_data$GroupCellType <- interaction(filtered_data$CellType, filtered_data$Group)
filtered_data <- set_group_celltype_order(filtered_data, selected_celltypes)

# Ensure CellType factor follows the same order
filtered_data$CellType <- factor(filtered_data$CellType, levels = selected_celltypes)

# Define color palette
mypalette <- c(
  "Epithelial" = "#A6CEE3",
  "Neuroendocrine" = "#FB9A99",
  "Fibroblast" = "#1F78B4",
  "Smooth muscle" = "#E31A1C",
  "Endothelial" = "#CAB2D6",
  "T & NK cells" = "darkred",
  "B cells" = "#33A02C",
  "Myeloblast" = "#FDBF6F",
  "Macrophage" = "darkgrey",
  "DC" = "#7570B3",
  "pDC" = "#D95F02",
  "Neutrophils-1" = "#6A3D9A",
  "Neutrophils-2" = "#B15928",
  "Eosinophils" = "#1B9E77",
  "Mast cells" = "#FF7F00",
  "Megakaryocyte" = "navy",
  "Erythroid cells" = "black"
)

# Prepare stat annotations for ggpubr
stat_annotations <- dunn_results %>%
  filter(CellType %in% selected_celltypes) %>%
  mutate(
    group1 = paste(CellType, group1, sep = "."),
    group2 = paste(CellType, group2, sep = "."),
    # Convert p-values to significance stars
    p.signif = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01 ~ "**",
      P.adj < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  )

# Calculate max values for each cell type to position brackets
max_values <- filtered_data %>%
  group_by(CellType) %>%
  summarise(max_prop = max(Proportion, na.rm = TRUE))

# Add y.position based on cell type
stat_annotations <- stat_annotations %>%
  left_join(max_values, by = "CellType") %>%
  group_by(CellType) %>%
  mutate(
    y.position = max_prop + 0.05 + (row_number() - 1) * 0.08
  ) %>%
  ungroup()

# Create base boxplot
p_dunn <- ggplot(filtered_data, aes(x = GroupCellType, y = Proportion, fill = CellType)) +
  geom_boxplot(outlier.shape = 21, outlier.color = "black") +
  geom_jitter(shape = 21, color = "black", size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = mypalette) +
  theme_bw(base_family = "Arial") +
  labs(x = NULL, y = "Proportion", fill = "Cell Types") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", family = "Arial", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold", family = "Arial"),
    axis.title = element_text(size = 14, face = "bold", family = "Arial"),
    legend.text = element_text(size = 10, face = "bold", family = "Arial"),
    legend.title = element_text(size = 12, face = "bold", family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 1.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = seq(3.5, 18, by = 3), linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_x_discrete(labels = rep(c("Control", "Asthma", "EGPA"), times = length(selected_celltypes))) +
  stat_pvalue_manual(
    stat_annotations,
    label = "p.signif",  # Use stars instead of p-values
    size = 3.8,
    label.size = 3.8,
    bracket.size = 0.5,
    tip.length = 0.01,
    step_increase = 0.01,
    hide.ns = TRUE  # Hide non-significant comparisons
  ) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, by = 0.2))

# Print the plot
print(p_dunn)

# Save plot
library(Cairo)
ggsave("./BALF_boxplot_with_Dunn_test.pdf",
       plot = p_dunn, width = 10, height = 6, device = cairo_pdf)
cat("\n=== Boxplot with Dunn's test saved! ===\n")

# ------------------------------------------------------------------------------
# Biopsy Analysis - Stacked Bar Plot
# ------------------------------------------------------------------------------
tissue <- "Biopsy"
data <- pdata[which(pdata$Var2 == tissue), ]
data <- data[which(!is.na(data$Proportion)), ]

# Calculate normalized proportions per group
data_summary <- data %>%
  group_by(Group, Var3) %>%
  summarise(Proportion = sum(Proportion, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Proportion = Proportion / sum(Proportion)) %>%
  ungroup()

# Recode group labels
data_summary$Group <- factor(data_summary$Group,
                             levels = c("E", "A", "C"),
                             labels = c("EGPA", "Asthma", "Control"))

# Create stacked bar plot
ggplot(data_summary, aes(x = Group, y = Proportion, fill = Var3)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_family = "Arial") +
  scale_fill_manual(values = mypalette) +
  labs(x = "", y = "Proportion", fill = "") +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", family = "Arial", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14, face = "bold", family = "Arial"),
    axis.title = element_text(size = 16, face = "bold", family = "Arial"),
    legend.text = element_text(size = 14, face = "bold", family = "Arial"),
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),
    panel.border = element_rect(color = "black", size = 1.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.ticks.length = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(ncol = 6, keywidth = 0.8, keyheight = 0.8)) +
  coord_flip()

# ------------------------------------------------------------------------------
# Biopsy Analysis - Box Plots with Statistics
# ------------------------------------------------------------------------------
tissue <- "Biopsy"
data <- pdata[which(pdata$Var2 == tissue), ]
data <- data[which(!is.na(data$Proportion)), ]

# Recode group labels
data <- data %>%
  mutate(Group = recode(Group,
                        "C" = "Control",
                        "A" = "Asthma",
                        "E" = "EGPA")) %>%
  dplyr::rename(CellType = Var3)

# Set factor levels
data$Group <- factor(data$Group, levels = c("Control", "Asthma", "EGPA"))

# ------------------------------------------------------------------------------
# Kruskal-Wallis Test with Effect Size (η²H)
# ------------------------------------------------------------------------------

# Function to calculate Kruskal-Wallis with effect size
perform_kw_test <- function(data) {
  # Get valid cell types (present in all 3 groups)
  valid_celltypes <- data %>%
    filter(Proportion > 0) %>%
    group_by(CellType) %>%
    summarise(n_groups = n_distinct(Group)) %>%
    filter(n_groups >= 3) %>%
    pull(CellType)
  
  # Filter data for valid cell types
  filtered_data <- data %>%
    filter(CellType %in% valid_celltypes)
  
  # Perform Kruskal-Wallis test for each cell type
  kw_results <- filtered_data %>%
    group_by(CellType) %>%
    kruskal_test(Proportion ~ Group) %>%
    mutate(
      # Calculate eta-squared H (effect size)
      eta_squared_H = statistic / (n - 1)
    )
  
  return(kw_results)
}

kw_results <- perform_kw_test(data)

# Print Kruskal-Wallis results
cat("\n=== Kruskal-Wallis Test Results ===\n")
print(kw_results)

# ------------------------------------------------------------------------------
# Dunn's Post-hoc Test with Bonferroni Adjustment
# ------------------------------------------------------------------------------

# Function to perform Dunn's test
perform_dunn_test <- function(data) {
  # Get valid cell types
  valid_celltypes <- data %>%
    filter(Proportion > 0) %>%
    group_by(CellType) %>%
    summarise(n_groups = n_distinct(Group)) %>%
    filter(n_groups >= 3) %>%
    pull(CellType)
  
  # Filter data for valid cell types
  filtered_data <- data %>%
    filter(CellType %in% valid_celltypes)
  
  # Perform Dunn's test for each cell type
  dunn_results_list <- list()
  
  for (cell_type in valid_celltypes) {
    cell_data <- filtered_data %>%
      filter(CellType == cell_type)
    
    # Perform Dunn's test using FSA package
    dunn_test <- dunnTest(Proportion ~ Group, 
                          data = cell_data, 
                          method = "bh")
    
    # Extract results
    dunn_df <- dunn_test$res %>%
      mutate(CellType = cell_type) %>%
      select(CellType, Comparison, Z, P.unadj, P.adj)
    
    dunn_results_list[[cell_type]] <- dunn_df
  }
  
  # Combine all results
  dunn_results <- bind_rows(dunn_results_list)
  
  # Separate comparison column into group1 and group2
  dunn_results <- dunn_results %>%
    separate(Comparison, into = c("group1", "group2"), sep = " - ", remove = FALSE)
  
  return(dunn_results)
}

dunn_results <- perform_dunn_test(data)

# Print Dunn's test results
cat("\n=== Dunn's Post-hoc Test Results (Bonferroni Adjusted) ===\n")
print(dunn_results)

# ------------------------------------------------------------------------------
# Create Comprehensive Results Table
# ------------------------------------------------------------------------------

# Merge KW and Dunn's results
# First convert kw_results to a regular data frame
kw_df <- as.data.frame(kw_results) %>%
  select(CellType, statistic, df, p, n, eta_squared_H) %>%
  dplyr::rename(KW_chi_squared = statistic,
                KW_df = df,
                KW_p_value = p,
                KW_n = n)

comprehensive_results <- dunn_results %>%
  left_join(kw_df, by = "CellType") %>%
  mutate(
    # Add significance levels
    significance = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01 ~ "**",
      P.adj < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    # Format p-values
    P.adj_formatted = case_when(
      P.adj < 0.001 ~ "< 0.001",
      TRUE ~ sprintf("%.4f", P.adj)
    )
  ) %>%
  arrange(CellType, group1, group2)

# Print comprehensive results
cat("\n=== Comprehensive Statistical Results ===\n")
print(comprehensive_results)

# ------------------------------------------------------------------------------
# Save Results
# ------------------------------------------------------------------------------

# Save Kruskal-Wallis results
write.csv(kw_results, 
          "./Biopsy_KW_test_results.csv", 
          row.names = FALSE)

# Save Dunn's test results
write.csv(dunn_results, 
          "./Biopsy_Dunn_test_results.csv", 
          row.names = FALSE)

# Save comprehensive results (for Supplementary Table 5)
write.csv(comprehensive_results, 
          "./Biopsy_Supplementary_Table5_statistics.csv", 
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Create Publication-Ready Table
# ------------------------------------------------------------------------------

# Format for publication
pub_table <- comprehensive_results %>%
  select(CellType, Comparison, Z, P.adj, P.adj_formatted, significance, 
         KW_chi_squared, KW_p_value, eta_squared_H) %>%
  mutate(
    # Round values for presentation
    Z = round(Z, 3),
    KW_chi_squared = round(KW_chi_squared, 3),
    eta_squared_H = round(eta_squared_H, 3),
    KW_p_value = case_when(
      KW_p_value < 0.001 ~ "< 0.001",
      TRUE ~ sprintf("%.4f", KW_p_value)
    )
  )

# Save publication-ready table
write.csv(pub_table, 
          "./Biopsy_Publication_Ready_Statistics.csv", 
          row.names = FALSE)

cat("\n=== Results saved successfully! ===\n")
cat("Files saved:\n")
cat("1. Biopsy_KW_test_results.csv\n")
cat("2. Biopsy_Dunn_test_results.csv\n")
cat("3. Biopsy_Supplementary_Table5_statistics.csv\n")
cat("4. Biopsy_Publication_Ready_Statistics.csv\n")

# ------------------------------------------------------------------------------
# Create Summary Report for Figure Legend
# ------------------------------------------------------------------------------

# Extract significant results for figure legend
significant_results <- comprehensive_results %>%
  filter(P.adj < 0.05) %>%
  select(CellType, group1, group2, P.adj_formatted, significance)

# Format for figure legend text
figure_legend_text <- significant_results %>%
  mutate(
    text = sprintf("%s: %s vs %s, P = %s", 
                   CellType, group1, group2, P.adj_formatted)
  ) %>%
  pull(text)

# Save figure legend text
writeLines(figure_legend_text, "./Biopsy_Figure_Legend_P_values.txt")

cat("\n=== Figure Legend P-values (P < 0.05) ===\n")
for (text in figure_legend_text) {
  cat(paste0("• ", text, "\n"))
}

# ------------------------------------------------------------------------------
# Create Enhanced Boxplot with Dunn's Test Results
# ------------------------------------------------------------------------------

# Filter for selected cell types (in the desired order)
selected_celltypes <- c("pDC", "Eosinophils", "DC", "Neutrophils-2", "Epithelial")

filtered_data <- data %>%
  filter(CellType %in% selected_celltypes)

# Function to set GroupCellType order (using custom order)
set_group_celltype_order <- function(data, cell_type_order) {
  groups <- c("Control", "Asthma", "EGPA")
  sorted_levels <- c()
  for (cell_type in cell_type_order) {
    for (group in groups) {
      sorted_levels <- c(sorted_levels, paste(cell_type, group, sep = "."))
    }
  }
  data$GroupCellType <- factor(data$GroupCellType, levels = sorted_levels)
  return(data)
}

# Create GroupCellType column and set order
filtered_data$GroupCellType <- interaction(filtered_data$CellType, filtered_data$Group)
filtered_data <- set_group_celltype_order(filtered_data, selected_celltypes)

# Ensure CellType factor follows the same order
filtered_data$CellType <- factor(filtered_data$CellType, levels = selected_celltypes)

# Define color palette
mypalette <- c(
  "Epithelial" = "#A6CEE3",
  "Neuroendocrine" = "#FB9A99",
  "Fibroblast" = "#1F78B4",
  "Smooth muscle" = "#E31A1C",
  "Endothelial" = "#CAB2D6",
  "T & NK cells" = "darkred",
  "B cells" = "#33A02C",
  "Myeloblast" = "#FDBF6F",
  "Macrophage" = "darkgrey",
  "DC" = "#7570B3",
  "pDC" = "#D95F02",
  "Neutrophils-1" = "#6A3D9A",
  "Neutrophils-2" = "#B15928",
  "Eosinophils" = "#1B9E77",
  "Mast cells" = "#FF7F00",
  "Megakaryocyte" = "navy",
  "Erythroid cells" = "black"
)

# Prepare stat annotations for ggpubr
stat_annotations <- dunn_results %>%
  filter(CellType %in% selected_celltypes) %>%
  mutate(
    group1 = paste(CellType, group1, sep = "."),
    group2 = paste(CellType, group2, sep = "."),
    # Convert p-values to significance stars
    p.signif = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01 ~ "**",
      P.adj < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  )

# Calculate max values for each cell type to position brackets
max_values <- filtered_data %>%
  group_by(CellType) %>%
  summarise(max_prop = max(Proportion, na.rm = TRUE))

# Add y.position based on cell type
stat_annotations <- stat_annotations %>%
  left_join(max_values, by = "CellType") %>%
  group_by(CellType) %>%
  mutate(
    y.position = max_prop + 0.05 + (row_number() - 1) * 0.08
  ) %>%
  ungroup()

# Create base boxplot
p_dunn <- ggplot(filtered_data, aes(x = GroupCellType, y = Proportion, fill = CellType)) +
  geom_boxplot(outlier.shape = 21, outlier.color = "black") +
  geom_jitter(shape = 21, color = "black", size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = mypalette) +
  theme_bw(base_family = "Arial") +
  labs(x = NULL, y = "Proportion", fill = "Cell Types") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", family = "Arial", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold", family = "Arial"),
    axis.title = element_text(size = 14, face = "bold", family = "Arial"),
    legend.text = element_text(size = 10, face = "bold", family = "Arial"),
    legend.title = element_text(size = 12, face = "bold", family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 1.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = seq(3.5, 14, by = 3), linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_x_discrete(labels = rep(c("Control", "Asthma", "EGPA"), times = length(selected_celltypes))) +
  stat_pvalue_manual(
    stat_annotations,
    label = "p.signif",  # Use stars instead of p-values
    size = 3.8,
    label.size = 3.8,
    bracket.size = 0.5,
    tip.length = 0.01,
    step_increase = 0.01,
    hide.ns = TRUE  # Hide non-significant comparisons
  ) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, by = 0.2))

# Print the plot
print(p_dunn)

# Save plot
library(Cairo)
ggsave("./Biopsy_boxplot_with_Dunn_test.pdf",
       plot = p_dunn, width = 10, height = 6, device = cairo_pdf)
cat("\n=== Boxplot with Dunn's test saved! ===\n")

# ==============================================================================
# End of Script
# ==============================================================================
