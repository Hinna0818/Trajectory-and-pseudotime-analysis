## change seurat object into sce object
# raw data from https://github.com/zhanghao-njmu/SCP/blob/b9b0eb7a7bf2c2c4b2262e73e09d7ebd515c7da0/data/pancreas_sub.rda
# script https://zhanghao-njmu.github.io/SCP/index.html

# load data
library(SingleCellExperiment)
rm(list = ls())
load("pancreas_sub.rda")
srt <- pancreas_sub

# extract counts_matrix and logcounts_matrix
counts_matrix <- as.matrix(srt@assays$RNA@counts)
logcounts_matrix <- as.matrix(srt@assays$RNA@data)

# extract cell_metadata
cell_metadata <- srt@meta.data

# extract gene_metadata
gene_metadata <- data.frame(gene_short_name = rownames(counts_matrix))
rownames(gene_metadata) <- rownames(counts_matrix)

# create sce object
sce <- SingleCellExperiment(
  assays = list(counts = counts_matrix, logcounts = logcounts_matrix),
  colData = cell_metadata,
  rowData = gene_metadata
)

# add reduction result
reducedDims(sce)$PCA <- as.matrix(Embeddings(srt, "PCA"))
reducedDims(sce)$UMAP <- as.matrix(Embeddings(srt, "UMAP"))

saveRDS(sce, file = "pancreas_sub_sce.rds")
