# dreamlet universe
library(dreamlet)
library(crumblr)

# data IO
library(SingleCellExperiment) 
library(zellkonverter)
library(tidyr)

# plotting
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(ggtree)
library(aplot)

# versions
sessionInfo()

prefix = 'example_run'

### Model Definition
# Define the mixed-effects model for differential expression
# - log(n_counts): Covariate to correct for sequencing depth (normalization)
# - (1|batch): Random effect to capture batch-to-batch variability
# - scale(age): Standardized age (fixed effect)
# - (1|sex) + (1|ancestry): Random effects for demographics
# - PMI + Braak: Post-mortem interval and disease stage (fixed effects of interest)
form = ~ log(n_counts) + (1|batch) + scale(age) + (1|sex) + (1|ancestry) + PMI + Braak + 1

### 1. Load Pseudobulk Object
# Reads the SingleCellExperiment object containing pseudobulk counts (aggregated)
pbObj <- readRDS(paste0(prefix,'.pbObj.rds'))

### 2. Processing & Normalization
# 'processAssays' performs filtering and normalization for each cell type
# - min.cells=5: Exclude genes present in <5 cells
# - min.count=5: Exclude genes with <5 counts
# - min.samples=4: Exclude genes present in <4 samples
# - min.prop=0.2: Exclude genes present in <20% of samples (optional)
# - Normalization: voom (or voomWithDreamWeights) to prep for linear modeling
res.proc = processAssays(pbObj, form, min.cells=5, min.count=5, min.samples=4, min.prop=0.2)
saveRDS(res.proc, paste0(prefix,'.res.proc.rds'))

### 3. Variance Partitioning
# Quantify the fraction of expression variation explained by each variable in 'form'
# This helps identify if Batch, PMI, or QC metrics are driving the data more than Biology (Braak)
vp.lst = fitVarPart(res.proc, form)
saveRDS(vp.lst, paste0(prefix,'.vp.lst.rds'))

# Quick check: How many genes passed filtering?
length(unique(vp.lst$gene))

# Aggregate variance fractions across all genes (mean)
vp.agg <- aggregate(. ~ gene, vp.lst[,-1], mean)
rownames(vp.agg) <- vp.agg$gene
vp.agg <- vp.agg[,-1]

# Plot: Visual summary of drivers of variation
# Look for biological variables explaining a significant portion of variance
options(repr.plot.width=5, repr.plot.height=5)
plotVarPart(sortCols(vp.agg), label.angle=45, ncol=1) + theme(aspect.ratio=1)

### 4. Differential Expression (Dreamlet)
# Fit the linear mixed model for each gene in each cell type
# Tests the effect of the variables in 'form' (specifically focusing on fixed effects like Braak)
res.dl.braak = dreamlet(res.proc, form)
saveRDS(res.dl.braak, paste0(prefix,'.res.dl.braak.rds'))
