# sc
import pegasus as pg
import scanpy as sc
import anndata as ad
import pge

# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# data
import numpy as np
import pandas as pd
from scipy import stats
from scipy import sparse
import h5py

### Analysis Parameters
# 'example_run' prefix used for all output files
prefix='example_run'
flavor='cell_ranger'     # HVF selection method (e.g., 'cell_ranger', 'seurat')
batch_key='Source'       # Column in .obs to use for batch correction/integration
n_top_genes=1000         # Number of highly variable features (HVGs) to select
npc=50                   # Number of Principal Components (PCs) to compute
K=100                    # Number of neighbors for kNN graph construction
n_neighbors=50           # Number of neighbors for UMAP/tSNE visualization
res=2.0                  # High resolution for granular clustering (Leiden algorithm)

### 1. Load Pre-processed Data
# Expects an AnnData object that has already undergone QC and normalization
data = pg.read_input(prefix+'_qc_norm.h5ad')

### 2. Feature Selection
# Remove genes flagged as 'non-robust' during QC (e.g., expressed in <0.05% of cells)
data._inplace_subset_var(data.var['robust'])

# Select Highly Variable Features (HVF) to drive the downstream analysis
# robust_protein_coding=True ensures we focus on biologically relevant coding genes, ignoring noise
pge.scanpy_hvf(data, flavor=flavor, batch_key=batch_key, n_top_genes=n_top_genes, robust_protein_coding=True)

### 3. Dimensionality Reduction & Integration
# Initial PCA on the selected HVFs
pg.pca(data, n_components=npc)

# Optional: Check the elbow plot to verify if npc=50 is appropriate (uncomment to view)
pg.elbowplot(data)

# Update npc if internal selection logic changed it
npc = data.uns["pca_ncomps"]

# Regress out technical covariates that might drive unwanted variation
# - n_counts: Sequencing depth
# - percent_mito: Cell stress/quality (mitochondrial content)
# - cycle_diff: Cell cycle phase effects (G2M vs S score)
pg.regress_out(data, attrs=['n_counts','percent_mito','cycle_diff'])

# Run Harmony to integrate data across batches defined by 'batch_key'
# Uses the regressed PCA as input to correct batch effects
pg.run_harmony(data, batch=batch_key, rep='pca_regressed', max_iter_harmony=20, n_comps=npc)

### 4. Graph Construction & Visualization
# Build kNN graph using the Harmony-corrected PCA space (cosine distance metric)
pg.neighbors(data, rep='pca_regressed_harmony', use_cache=False, dist='cosine', K=K, n_comps=npc)

# Compute UMAP and tSNE embeddings for visualization
pg.umap(data, rep='pca_regressed_harmony', n_neighbors=n_neighbors, rep_ncomps=npc)
pg.tsne(data, rep='pca_regressed_harmony', rep_ncomps=npc)

### 5. Clustering (Leiden)
# Perform community detection using the Leiden algorithm.
# Resolution 2.0 yields a finer-grained clustering.
# Resulting labels are stored in .obs['leiden_labels_res20']
pg.leiden(data, rep='pca_regressed_harmony', resolution=res, class_label='leiden_labels_res20')

# Differential Expression (DE) Analysis
# Identify marker genes for each cluster against the rest
pg.de_analysis(data, cluster='leiden_labels_res20')

# Extract markers and save DE results to Excel for manual review
marker_dict_res20 = pg.markers(data)
pg.write_results_to_excel(marker_dict_res20, prefix+"_pass2_res20.de.xlsx")

# Auto-Annotation
# Infer cell types using a predefined marker list ('human_brain')
celltype_dict_res20 = pg.infer_cell_types(data, markers='human_brain')
cluster_names_res20 = pg.infer_cluster_names(celltype_dict_res20)

# Add the inferred annotations to the data object
pg.annotate(data, name='anno_res20', based_on='leiden_labels_res20', anno_dict=cluster_names_res20)

### 6. Visualization
# Generate scatter plots colored by Leiden clusters and inferred annotations
# legend_loc='on data' places labels directly on the plot
pg.scatter(data, attrs=['leiden_labels_res20','anno_res20'], basis='umap', legend_loc='on data', wspace=0.1)
