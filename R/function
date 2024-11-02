#' wrapper function for running monocle3 on default parameter
#' 
#' https://cole-trapnell-lab.github.io/monocle3/
#' 
#' @param sce a singlecellExperiment object.
#' @param cluster_coloumn the coloumn when clustering in the entered sce object.
#' @param assay_name the name of assay.
#' @param reduction_method the reduction used, using UMAP as default.
#' @param batch whether to remove batch effect, user can choose it personally.
#' @param seed The random seed to use for reproducibility. Defaults to 2024-11-1.
#' 
#' @importFrom SummarizedExperiment assay.
#' @importFrom SummarizedExperiment assayNames.

#' @export


runMonocle <- function(
    sce,
    assay_name = "logcounts", ## use logcounts as default
    reduction_method = "UMAP", ## use UMAP as default reduction method
    cluster_coloumn = "clusters", ## use clusters as default when clustering
    use_preprocess = FALSE, ## not preprocess the data as default, but can use preprocessed data in the previous pipeline
    batch = FALSE, ## not remove batch effects as default, but can do it personally
    num_dims = 50, ## set dims in PCA as 50 as default
    resolution = 1e-5, ## set resolution in cell-clustering
    seed = 2024-11-1
    ){
  
  set.seed(seed)
  
  ## 0. ensure the data is correctly set
  if(!"logcounts" %in% assayNames(sce)){
    print("logcounts assay does not exist in the object")
    sce <- scuttle::logNormCounts(sce)  ## use sclet function instead
  }
  
  ## 1. set data
  counts_matrix <- SummarizedExperiment::assay(sce, assay_name)
  gene_metadata <- data.frame(gene_short_name = rownames(counts_matrix))
  cell_metadata <- as.data.frame(colData(sce))
  rownames(gene_metadata) <- rownames(counts_matrix)
  
  ## create cellDataset
  cds <- monocle3::new_cell_data_set(counts_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_metadata)
  
  
  ## 2. preprocess the data
  if(use_preprocess){
    # Normalize and preprocess the data using preprocess_cds in Monocle3
    cds <- monocle3::preprocess_cds(cds, num_dim = num_dims)
    
    # create elbow plot, then choose the optimal number of dims
    print("Generating elbow plot to choose optimal number of dimensions...")
    print(monocle3::plot_pc_variance_explained(cds))
    
    # let user enter new dims number
    message("Please inspect the elbow plot to choose the optimal number of dimensions.")
    new_num_dims <- as.numeric(readline(prompt = "Enter the number of dimensions to use: "))
    
    # overlapped the former input
    if (!is.na(new_num_dims) && new_num_dims != num_dims) {
      num_dims <- new_num_dims
      cds <- monocle3::preprocess_cds(cds, num_dim = num_dims)
    }
    
    # Reduce dimensionality and visualize the cells
    cds <- monocle3::reduce_dimension(cds, reduction_method = reduction)
    
    # move batch effect personally
    while(batch && "plate" %in% colnames(colData(cds))){
      cds <- monocle3::align_cds(cds, num_dim = num_dims, alignment_group = "plate")
      cds <- reduce_dimension(cds)
    }
    
    # Group cells into cluster
    cds <- monocle3::cluster_cells(cds, resolution = resolution)
  }
  
  else{
    ## user can use their own pre-processed data instead of monocle3 default
    if(!is.null(SingleCellExperiment::reducedDims(sce))){
      reducedDims(cds)$UMAP <- reducedDims(sce)[[reduction_method]] ## using UMAP as default
    }
    else{
      stop("The object you entered does not contain reduced dimensions")
    }
    
    if(cluster_coloumn %in% colnames(colData(sce))){
      cds@clusters$UMAP$clusters <- factor(cell_metadata[[cluster_coloumn]])
    }
    else{
      stop("The object you entered does not contain specified cluster information")
    }
    
  }
  ## 3. trajectory analysis
  cds <- monocle3::learn_graph(cds)
  
  ## 4. pseudotime analysis
  cds <- monocle3::order_cells(cds)
  
  return(cds)
}
