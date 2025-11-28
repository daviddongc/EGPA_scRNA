#!/usr/bin/env python
# ==============================================================================
# EGPA scRNA-seq Analysis: Figure S4
# ==============================================================================
# Description: This script generates Supplementary Figure S4
# for <Airway immune profiles and therapeutic implications of IGF1 in eosinophilic
# granulomatosis with polyangiitis> publication
# ==============================================================================

import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy as scipy
import matplotlib.pyplot as plt
import gseapy as gp
import omicverse as ov
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

# Set matplotlib parameters
mpl.rcParams['pdf.fonttype'] = 42
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, color_map='Blues')

plt.rcParams.update({
    "axes.grid": False,
    "axes.facecolor": "white"
})

# ==============================================================================
# FIGURE S4N: GSEA Analysis for Eosinophil Subsets
# ==============================================================================
Eos_4_EA = sc.read('./data/FigS4/S4_n_EOS_4_EA.h5ad')
Eos_4_EA.obs['group'] = pd.Categorical(Eos_4_EA.obs['group'], categories=["E", "A"], ordered=True)
indices = Eos_4_EA.obs.sort_values(['leiden', 'group']).index
Eos_4_EA = Eos_4_EA[indices, :]

res = gp.gsea(data=Eos_4_EA.to_df().T,
              gene_sets="MSigDB_Hallmark_2020",
              cls=Eos_4_EA.obs['group'],
              permutation_num=1000,
              permutation_type='phenotype',
              outdir=None,
              method='s2n',
              threads=16)

terms = res.res2d.Term
axs = res.plot(terms[:5], show_ranking=True, legend_kws={'loc': (1.05, -0.55)})
plt.savefig("./data/FigS4/S4_n.pdf", dpi=300, format='pdf', bbox_inches="tight")

# ==============================================================================
# FIGURE S4M: Enrichment Analysis for Eosinophil-APOC1
# ==============================================================================
Eos_real = sc.read('./data/FigS4/S4_m_Eos_real.h5ad')
sc.tl.rank_genes_groups(Eos_real, groupby='subset_3', method='wilcoxon', pts=True)
result = Eos_real.uns['rank_genes_groups']
groups = result['names'].dtype.names
degs = pd.DataFrame(
    {group + '_' + key: result[key][group]
     for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
degs_sig = degs[degs['Eosinophil-APOC1_pvals_adj'] < 0.05]
degs_up = degs_sig[degs_sig['Eosinophil-APOC1_logfoldchanges'] > 0.5]
degs_dw = degs_sig[degs_sig['Eosinophil-APOC1_logfoldchanges'] < 0]
enr_up = gp.enrichr(degs_up['Eosinophil-APOC1_names'],
                    gene_sets='MSigDB_Hallmark_2020',
                    outdir=None)
colors = ['#ff0000', '#ff5500', '#ffaa00', '#fffe00']
cmap = LinearSegmentedColormap.from_list("autumn", colors[::-1])
gp.dotplot(enr_up.res2d, column='Adjusted P-value', figsize=(5, 7), title="Up", cmap=cmap, top_term=5)
plt.savefig("./data/FigS4/S4_m.pdf", dpi=300, bbox_inches='tight')

# ==============================================================================
# FIGURE S4L: Dotplot of Eosinophil Markers
# ==============================================================================
cell_dict = {'0': 'Eosinophil-SYNE2',
             '1': 'Eosinophil-CHI3L1',
             '3': 'Eosinophil-IL1RL1',
             '4': 'Eosinophil-APOC1'}
Eos_real.obs['subset_3'] = Eos_real.obs['leiden'].map(cell_dict)
eos_markers = [
    "SYNE2", "ROCK1", "ITGB2", "TLN1", "SIGLEC10",
    "CHI3L1", "CCL4L2", "ACP5", "ALOX5AP", "CCL18",
    "CCL4", "AREG", "CD69", "BCL2A1", "ALOX5",
    "APOC1", "TXNIP", "SELL",
]
desired_order = ['Eosinophil-SYNE2',
                 'Eosinophil-CHI3L1',
                 'Eosinophil-IL1RL1',
                 'Eosinophil-APOC1']
Eos_real.obs['subset_3'] = pd.Categorical(Eos_real.obs['subset_3'], categories=desired_order, ordered=True)
sc.pl.dotplot(Eos_real, var_names=eos_markers, groupby='subset_3', standard_scale='var', cmap='RdBu_r',
              categories_order=desired_order, dot_max=0.5, var_group_rotation=45,
              save='./data/FigS4/Eos_4_marker_dotplot.pdf')

# ==============================================================================
# FIGURE S4I: UMAP of Eosinophil Subsets
# ==============================================================================
custom_palette = sns.color_palette("Dark2", 4).as_hex()
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Eos_real,
                basis='X_umap',
                color='subset_3',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_I.pdf', dpi=600, bbox_inches='tight')

# ==============================================================================
# FIGURE S4J: UMAP by Group
# ==============================================================================
# Control
Eos_real_C = Eos_real[Eos_real.obs['group'].isin(['C']), :]
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Eos_real_C,
                basis='X_umap',
                color='subset_3',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_J1.pdf', dpi=600, bbox_inches='tight')

# Asthma
Eos_real_A = Eos_real[Eos_real.obs['group'].isin(['A']), :]
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Eos_real_A,
                basis='X_umap',
                color='subset_3',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_J2.pdf', dpi=600, bbox_inches='tight')

# EGPA
Eos_real_E = Eos_real[Eos_real.obs['group'].isin(['E']), :]
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Eos_real_E,
                basis='X_umap',
                color='subset_3',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_J3.pdf', dpi=600, bbox_inches='tight')

# ==============================================================================
# FIGURE S4K: Cell Proportion Analysis
# ==============================================================================
fig, ax = plt.subplots(figsize=(2, 7))
ov.pl.cellproportion(adata=Eos_real, celltype_clusters='subset_3',
                     groupby='group', legend=True, ax=ax)
plt.savefig('./data/FigS4/S4_K.pdf', dpi=300, bbox_inches='tight')

# ==============================================================================
# FIGURE S4O: Interferon Gene Expression
# ==============================================================================
colors = ["#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#E0F3F8", "#ABD9E9"]
cmap = LinearSegmentedColormap.from_list("RdYlBu_rev", colors[::-1])
interferon_related_genes = ["IFIT1", "IFIT2", "ISG15", "ISG15", "STAT1", "STAT2", "MX1", "RSAD2"]

# Load UMAP coordinates
umap_coords = pd.read_csv('./data/FigS4/S4_umap_coordinates.csv', index_col=0)
Eos_real.obsm['X_umap'] = umap_coords.loc[Eos_real.obs_names, ['UMAP_1', 'UMAP_2']].values

sc.pl.umap(Eos_real, color=interferon_related_genes, cmap=cmap, size=8,
           save='./data/FigS4/S4_O.pdf')

# ==============================================================================
# FIGURE S4H: Mast and Basophil Analysis
# ==============================================================================
Ma_Ba_single_real = sc.read('./data/FigS4/S4_H_Ma_Ba_single_real_sub6.h5ad')
sc.tl.rank_genes_groups(Ma_Ba_single_real, groupby='leiden_sub6', method='wilcoxon', pts=True)
ba = Ma_Ba_single_real[Ma_Ba_single_real.obs['leiden_sub6'].isin(['3,0']), :]
genes = [
    "IFITM3", "IFITM1", "IFIT1", "IFIT2", "RSAD2",
    "STAT1", "STAT2", "TYK2", "OAS1", "OAS2",
    "OAS3", "MX1", "MX2", "ISG15",
    "IFI6", "IFI27",
    "BCL2A1", "SOCS3",
    "TRIM21",
]
desired_order = ['C', 'A', 'E']
sc.pl.matrixplot(ba, var_names=genes, groupby='group',
                 standard_scale='var', cmap='RdBu_r', categories_order=desired_order,
                 swap_axes=True,
                 var_group_rotation=45, save='./data/FigS4/S4_H.pdf')

# ==============================================================================
# FIGURE S4G: GSEA for Macrophage Subset
# ==============================================================================
ma1 = Ma_Ba_single_real[Ma_Ba_single_real.obs['leiden_sub6'].isin(['1']), :]
ma1_EA = ma1[ma1.obs['group'].isin(['E', 'A']), :]
ma1_EA.obs['group'] = pd.Categorical(ma1_EA.obs['group'], categories=["E", "A"], ordered=True)
indices = ma1_EA.obs.sort_values(['leiden_sub6', 'group']).index
ma1_EA = ma1_EA[indices, :]
res = gp.gsea(data=ma1_EA.to_df().T,
              gene_sets="MSigDB_Hallmark_2020",
              cls=ma1_EA.obs['group'],
              permutation_num=1000,
              permutation_type='phenotype',
              outdir=None,
              method='s2n',
              threads=16)
terms = res.res2d.Term
axs = res.plot(terms[:5], show_ranking=True, legend_kws={'loc': (1.05, -0.55)})
plt.savefig("./data/FigS4/S4_g.pdf", dpi=300, format='pdf', bbox_inches="tight")

# ==============================================================================
# FIGURE S4F: Enrichment for Subset 1
# ==============================================================================
result = Ma_Ba_single_real.uns['rank_genes_groups']
groups = result['names'].dtype.names
degs = pd.DataFrame(
    {group + '_' + key: result[key][group]
     for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
degs_sig = degs[degs['1_pvals_adj'] < 0.05]
degs_up = degs_sig[degs_sig['1_logfoldchanges'] > 1]
degs_dw = degs_sig[degs_sig['1_logfoldchanges'] < 0]
enr_up = gp.enrichr(degs_up['1_names'],
                    gene_sets='MSigDB_Hallmark_2020',
                    outdir=None)
gp.dotplot(enr_up.res2d, column='Adjusted P-value', figsize=(5, 7), title="Up", cmap=plt.cm.autumn_r,
           top_term=5, size=20)
plt.savefig("./data/FigS4/S4_f.pdf", dpi=300, bbox_inches='tight')

# ==============================================================================
# FIGURE S4E: Mast cell-FOS Enrichment
# ==============================================================================
result = Ma_Ba_single_real.uns['rank_genes_groups']
groups = result['names'].dtype.names
degs = pd.DataFrame(
    {group + '_' + key: result[key][group]
     for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})

degs_sig = degs[degs['Mast cell-FOS_pvals_adj'] < 0.05]
degs_up = degs_sig[degs_sig['Mast cell-FOS_logfoldchanges'] > 0]
degs_dw = degs_sig[degs_sig['Mast cell-FOS_logfoldchanges'] < 0]
enr_up = gp.enrichr(degs_up['Mast cell-FOS_names'],
                    gene_sets='MSigDB_Hallmark_2020',
                    outdir=None)
gp.dotplot(enr_up.res2d, column='Adjusted P-value', figsize=(5, 8), title="Mast cell-FOS DEGS Up",
           cmap=plt.cm.autumn_r, top_term=5, size=12)
plt.savefig("./data/FigS4/S4_e.pdf", dpi=300, bbox_inches='tight')

# ==============================================================================
# FIGURE S4D: Dotplot of Mast and Basophil Markers
# ==============================================================================
desired_order = ['Mast cell-TPSAB1',
                 'Mast cell-SCN7A',
                 'Mast cell-FOSB',
                 'Mast cell-Cycling',
                 'Basophil-IL3RA']
Ma_Ba_single_real.obs['subset_4'] = pd.Categorical(Ma_Ba_single_real.obs['subset_4'], categories=desired_order,
                                                    ordered=True)
ma_ba_markers = [
    "TPSAB1", "TPSB2", "TPSD1", "PI3", "PLAT", "SCIN", "SNHG7",
    "SCN7A", "SLC12A8", "SH2D2A", "MAG", "RGS16",
    "FOSB", "JUN", "EGR1", "RGS1", "CD69", "JUNB",
    "BIRC5", "TOP2A", "TYMS", "MKI67", "UBE2C",
    "IL3RA", "CCR3", "HDC", "SELL", "MS4A3", "CSF2RB",
]
sc.pl.dotplot(Ma_Ba_single_real, var_names=ma_ba_markers, groupby='subset_4', standard_scale='var', cmap='RdBu_r',
              categories_order=desired_order, dot_max=0.5, var_group_rotation=45,
              save='./data/FigS4/S4_d.pdf')

# ==============================================================================
# FIGURE S4A: UMAP of All Mast and Basophil
# ==============================================================================
custom_palette = ['#1b9e77', '#d95f02', '#e7298a', '#66a61e', '#7570b3']
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Ma_Ba_single_real,
                basis='X_umap',
                color='subset_4',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_a.pdf', dpi=600, bbox_inches='tight')

Ma_Ba_single_real_C = Ma_Ba_single_real[Ma_Ba_single_real.obs['group'].isin(['C']), :]
Ma_Ba_single_real_A = Ma_Ba_single_real[Ma_Ba_single_real.obs['group'].isin(['A']), :]
Ma_Ba_single_real_E = Ma_Ba_single_real[Ma_Ba_single_real.obs['group'].isin(['E']), :]

# ==============================================================================
# FIGURE S4B: UMAP by Group for Mast and Basophil
# ==============================================================================
# Control
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Ma_Ba_single_real_C,
                basis='X_umap',
                color='subset_4',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_b1.pdf', dpi=600, bbox_inches='tight')

# Asthma
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Ma_Ba_single_real_A,
                basis='X_umap',
                color='subset_4',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_b2.pdf', dpi=600, bbox_inches='tight')

# EGPA
fig, ax = plt.subplots(figsize=(4, 4))
ov.pl.embedding(Ma_Ba_single_real_E,
                basis='X_umap',
                color='subset_4',
                frameon='small',
                show=False,
                palette=custom_palette,
                size=20,
                ax=ax)
plt.title('Cell Type', fontsize=15)
plt.savefig('./data/FigS4/S4_b3.pdf', dpi=600, bbox_inches='tight')

# ==============================================================================
# FIGURE S4C: Cell Proportion for Mast and Basophil
# ==============================================================================
fig, ax = plt.subplots(figsize=(2, 7))
ov.pl.cellproportion(adata=Ma_Ba_single_real, celltype_clusters='subset_4',
                     groupby='group', legend=True, ax=ax)
plt.savefig('./data/FigS4/S4_c.pdf', dpi=300, bbox_inches='tight')

# ==============================================================================
# End of Script
# ==============================================================================
