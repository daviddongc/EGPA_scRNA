# Airway Immune Profiles and Therapeutic Implications of IGF1 in Eosinophilic Granulomatosis with Polyangiitis

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-%3E%3D3.8-blue)](https://www.python.org/)

## Overview

This repository contains all analysis code used to generate figures for the manuscript:

**"Airway immune profiles and therapeutic implications of IGF1 in eosinophilic granulomatosis with polyangiitis"**

The code includes comprehensive workflows for single-cell RNA sequencing (scRNA-seq) analysis, immune cell profiling, differential gene expression analysis, pathway enrichment, and therapeutic target identification in eosinophilic granulomatosis with polyangiitis (EGPA).

---

## Repository Structure

### Main Figures

```
Fig1.R          # Figure 1C-D: UMAP visualization and cell type composition
                # - Single-cell landscape of airway immune cells
                # - Cell type proportions in BALF and Biopsy samples
                # - Statistical comparisons between EGPA, Asthma, and Control groups

Fig2.R          # Figure 2A-F: Granulocyte subset analysis
                # - UMAP clustering of granulocyte populations
                # - Cellular source proportions across samples
                # - Feature plots for inflammatory and interferon-response genes
                # - Differential gene expression in circulating vs tissue eosinophils

Fig3.R          # Figure 3: B cell and Plasma cell analysis
                # - RNA velocity analysis (3A)
                # - Marker gene expression feature plots (3B)
                # - Cellular source proportions for B/Plasma subtypes (3C)
                # - Fc receptor expression in granulocytes (3F)

Fig4.R          # Figure 4: IGF1+ macrophage identification and characterization
                # - UMAP split by treatment status (4A)
                # - Volcano plot for IGF1 macrophage markers (4B)
                # - GO/Reactome pathway enrichment (4C)
                # - IGF1, HP, RETN expression feature plots (4D)
                # - MAGIC-imputed violin plots for pDC and macrophage (4E)
                # - Epithelial cell UMAP by treatment (4G)
                # - Relative enrichment heatmap (4H)
                # - Hormone receptor expression in epithelial cells (4I)
                # - MAGIC-imputed violin plots for Goblet and Ionocyte cells (4J)
```

### Supplementary Figures

```
FigS2.R         # Supplementary Figure 2E-G: Quality control and batch assessment
                # - UMAP colored by source, group, and subject
                # - Gene expression correlation analysis
                # - Heatmap of differentially expressed genes

FigS3.R         # Supplementary Figure 3: Extended granulocyte analysis
                # - Marker gene heatmap (S3A)
                # - SCENIC transcription factor activity heatmap (S3B)
                # - Granule gene expression violin plots (S3C)
                # - Ligand-receptor dotplots (S3D): CSF, Chemokine, Interleukin

FigS4.py        # Supplementary Figure 4: Mast/Basophil and Eosinophil subsets (Python)
                # - Mast cell and Basophil UMAP and markers (S4_a-e)
                # - Eosinophil subset characterization (S4_I-O)
                # - GSEA enrichment analysis
                # - Interferon-response gene visualization

FigS5.R         # Supplementary Figure 5: Circulatory vs Tissue eosinophils
                # - Differential expression heatmap
                # - GO enrichment bar plots (Metascape)

FigS7.R         # Supplementary Figure 7: T/NK and B cell immune checkpoint analysis
                # - T/NK cell UMAP and proportions (S7_a-b)
                # - Co-stimulatory/inhibitory molecule expression
                # - Complement pathway ligand expression
                # - B cell ligand-receptor dotplots

FigS9.R         # Supplementary Figure 9: Treatment response in macrophage/DC populations
                # - After-treatment UMAP (S9B)
                # - MAGIC-imputed violin plots for Neutrophils-2 and Eosinophils (S9D)
                # - Cell proportion bar plots by treatment and disease status (S9E)
                # - Marker gene heatmap (S9F)
                # - Relative enrichment (Ro/e) heatmap (S9G)

FigS10.R        # Supplementary Figure 10: Epithelial and Fibroblast analysis
                # - Epithelial cell type proportions (S10_a)
                # - Fibroblast UMAP (S10_c)
                # - IGF1 and CCL11 feature plots in fibroblasts (S10_d)
                # - Adventitial fibroblast DEG analysis with MAST (S10_e)
                # - MAGIC-imputed violin plots with statistical testing (S10_f)

FigS11.R        # Supplementary Figure 11: IGF1 signaling receptor expression
                # - IGF1R, IL17RA, IL17RB expression across cell types
                # - Feature plots in combined dataset, granulocytes, and T/NK cells
```


---

## Input Data

**Note:** Due to patient privacy protection, raw sequencing data and clinical information are available from the corresponding author upon reasonable request and with appropriate institutional review board (IRB) approval.

### Required Input Data Formats:
- **Seurat objects** (`.RDS` files): Preprocessed single-cell data
- **h5ad files**: AnnData objects for Python-based analysis
- **Metadata** (`.txt`, `.csv`): Sample information, clinical annotations
- **Expression matrices**: Gene expression count tables

### Data Organization:
```
data/
├── combined.RDS                  # All cells, integrated dataset
├── granu.RDS                     # Granulocyte subset
├── macro.RDS                     # Macrophage subset (baseline)
├── macroTreat.RDS               # Macrophage subset (post-treatment)
├── Bsub.RDS                     # B cell and Plasma cell subset
├── TNK.RDS                      # T & NK cell subset
├── fibro.RDS                    # Fibroblast subset
├── S4_*.h5ad                    # Python-processed data for S4
└── metadata.txt                  # Sample metadata
```

---

## Output

### Figure Formats
Scripts generate publication-ready figures in multiple formats:
- **PDF** (vector format, 300 dpi) - recommended for manuscripts
- **PNG** (raster format, 80-600 dpi) - for presentations/previews

### Output Files
- **Figures**: Saved in working directory or `figures/` subdirectories
- **Tables**:
  - Differential expression results (`.xls`, `.csv`)
  - Statistical test results
  - Enrichment analysis results
- **R objects**: Intermediate processed data (`.RDS`, `.rdata`)

---

## Citation

If you use this code or find our work helpful, please cite:

```
[Author List]. Airway immune profiles and therapeutic implications of IGF1
in eosinophilic granulomatosis with polyangiitis.
[Journal Name]. [Year]; [Volume]([Issue]):[Pages].
DOI: [DOI number]
```

---

## Contact

For questions, issues, or collaboration opportunities:

- **GitHub Issues**: [Open an issue](../../issues)
- **Email**: [Corresponding author email]
- **Institution**: [Your institution name]

---

## Acknowledgments

We thank all patients who participated in this study and the clinical staff who facilitated sample collection. We acknowledge the use of the following key computational resources:

- **Seurat** (Butler et al., 2018; Hao et al., 2021)
- **Scanpy** (Wolf et al., 2018)
- **Harmony** (Korsunsky et al., 2019)
- **MAST** (Finak et al., 2015)
- **MAGIC** (van Dijk et al., 2018)
- **SCENIC** (Aibar et al., 2017)
- **Velocyto** (La Manno et al., 2018)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---


**Last updated:** 2025-11-26
**Code version:** 1.0.0
