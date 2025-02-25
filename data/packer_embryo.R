## preprocessing data
##
## This script performs preprocessing on raw SingleCellExperiment data using sclet function for downstream analysis. 
##
## Input:
## - Raw gene expression data (exp) in the form of a matrix
## - Cell metadata (cell) in the form of a DataFrame
## - Gene annotation (gene) in the form of a DataFrame
## The original data is from https://github.com/cole-trapnell-lab/monocle3/blob/master/examples/c_elegans_embryo.R
##
## Output:
## - Processed SCE object saved as embryo_sce.rds


## load raw data
rm(list = ls())
exp <- readRDS("./packer_embryo_expression.rds")
cell <- readRDS("packer_embryo_colData.rds")
gene <- readRDS("packer_embryo_rowData.rds")

## load packages
library(SingleCellExperiment)
library(yulab.utils)
library(ggplot2)
source("./seurat.R")
source("./subset.R")

## create raw SCE object
embryo_sce <- SingleCellExperiment(assay = list(counts = exp), colData = cell, rowData = gene)
rownames(embryo_sce@assays@data$counts) <- gene[['gene_short_name']]
rownames(embryo_sce) <- gene[['gene_short_name']]

# QC
embryo_sce[["percent.mt"]] <- PercentageFeatureSet(embryo_sce, "^MT-")
subset_feature(embryo_sce, 10)
embryo_sce <- subset_feature(embryo_sce, mincell = 3, peek=FALSE)

# Normalization
embryo_sce <- NormalizeData(embryo_sce)

# correct batch effect using "batch" factor in ColData(embryo_sce)
embryo_sce_correct <- batchelor::quickCorrect(embryo_sce, batch = embryo_sce$batch)
corrected_logcounts <- assay(embryo_sce_correct$corrected, "reconstructed")
corrected_logcounts_dgc <- as(corrected_logcounts, "dgCMatrix")
assays(embryo_sce)$logcounts <- corrected_logcounts_dgc

# find variable features
embryo_sce <- FindVariableFeatures(embryo_sce)

# dimension reduction
embryo_sce <- ScaleData(embryo_sce)
embryo_sce <- scater::runPCA(embryo_sce, subset_row = VariableFeatures(embryo_sce), exprs_values = "scaled")
ElbowPlot(embryo_sce)

# clustering
set.seed(1234)
embryo_sce <- FindNeighbors(embryo_sce, dims = 1:20)
embryo_sce <- FindClusters(embryo_sce)

# UMAP
embryo_sce <- RunUMAP(embryo_sce, 1:20)
ggsc::sc_dim(embryo_sce, reduction = "UMAP") + ggsc::sc_dim_geom_label()

saveRDS(embryo_sce, "./embryo_sce.rds")
#usethis::use_data(embryo_sce, overwrite = TRUE)
