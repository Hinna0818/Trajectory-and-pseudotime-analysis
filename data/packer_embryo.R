## preprocessing the example data and convert it into SingleCellExperiment object
## The original data is from 

## load data
expression_matrix <- readRDS("packer_embryo_expression.rds")
cell_metadata <- readRDS("packer_embryo_colData.rds")
gene_annotation <- readRDS("packer_embryo_rowData.rds")

## create cds object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## preprocess data
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
cds <- reduce_dimension(cds)


## convert them into SCE object
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = expression_matrix),
                            colData = cell_metadata, 
                            rowData = gene_annotation)

reducedDims(sce)[["UMAP"]] <- reducedDim(cds, "UMAP")
rownames(sce@assays@data$counts) <- gene_annotation[['gene_short_name']]
rownames(sce) <- gene_annotation[['gene_short_name']]

saveRDS(sce, "packer_embryo.rds")
