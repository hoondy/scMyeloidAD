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
csv_manifest = '/path/to/manifest.csv'  # Path to the sample metadata/manifest file
list_attributes = ['Donor','Source','Type','Inhibitors','Sequencing'] # Columns from manifest to include in AnnData.obs
genome = 'GRCh38'                       # Reference genome (e.g., 'GRCh38', 'mm10')

# Quality Control Thresholds
n_genes_lower = 500       # Remove cells with < 500 genes detected
n_counts_lower = 1000     # Remove cells with < 1000 total counts (UMIs)
n_genes_upper = 8000      # Remove potential doublets (> 8000 genes)
n_counts_upper = 40000    # Remove potential doublets (> 40000 counts)
percent_cells = 0.05      # Gene must be detected in > 0.05% of cells to be considered robust
percent_mito = 20         # Remove cells with > 20% mitochondrial reads (dead/dying cells)

### 1. Data Aggregation
# Load raw count matrices from individual samples listed in the manifest
# - default_ref: genome reference used for mapping
# - append_sample_name: ensures unique cell barcodes across samples
data = pg.aggregate_matrices(csv_file=csv_manifest, attributes=list_attributes, default_ref=genome, append_sample_name=True)

### 2. Identify Robust Genes
# Flag genes expressed in at least 'percent_cells' of the dataset
# This helps filter out noisy, lowly expressed genes before downstream analysis
pg.identify_robust_genes(data, percent_cells=percent_cells)

### 3. Feature Selection (Filtering)
# Remove features that are not robust OR are not protein-coding (if applicable)
# Keeps the dataset clean and focused on informative genes
data._inplace_subset_var(data.var['robust'])

### 4. Quality Control (QC)
# Calculate QC metrics (n_genes, n_counts, percent_mito)
# - mito_prefix='MT-': prefix for mitochondrial genes in human (use 'mt-' for mouse)
pg.qc_metrics(data, min_genes=n_genes_lower, max_genes=n_genes_upper, min_umis=n_counts_lower, max_umis=n_counts_upper, mito_prefix='MT-', percent_mito=percent_mito)

# Apply the QC filters to remove low-quality cells and doublets
pg.filter_data(data)

### 5. Cleanup
# Remove unused categories from categorical columns (e.g., if a donor was filtered out entirely)
data = pge.clean_unused_categories(data)

### 6. Normalization
# Log-normalize the count data (Natural Log(CPM + 1))
# CPM = Counts Per Million
pg.log_norm(data)
