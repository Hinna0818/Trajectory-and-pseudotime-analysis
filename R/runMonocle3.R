#' Wrapper function for running monocle3 with default parameters.
#' 
#' @param sce SingleCellExperiment object.
#' @param assay_name The name of assay. Defaults to "logcounts".
#' @param use_auto_root Choose roots automatically in pseudotime analysis. Defaults to FALSE.
#' @param time_bin Character. The time bin to use for automatic root selection, e.g., "130-170".
#' @param reduction_method Dimensionality reduction method to use, default is "UMAP".
#' @param resolution The resolution parameter for clustering. Defaults to 1e-5.
#' @param cluster_method The clustering method to use. Defaults to "louvain".
#' @param k The number of nearest neighbors to consider for clustering. Defaults to 20.
#' @param num_iter The number of iterations for clustering. Defaults to 2.
#' @param seed Random seed for reproducibility.
#' 
#' @importFrom SummarizedExperiment assay assayNames
#' @export

runMonocle <- function(
    sce,
    assay_name = "logcounts",
    use_auto_root = FALSE,
    time_bin = "130-170",
    reduction_method = "UMAP",
    resolution = 1e-5,
    cluster_method = "louvain",
    k = 20,
    seed = 2024-11-01){

  set.seed(seed)
  
  reduction_method <- match.arg(reduction_method, c("UMAP", "tSNE", "PCA"))
  cluster_method <- match.arg(cluster_method, c("leiden", "louvain"))
  
  ## check the data before setting a cds object
  # ensure the logcounts assay exists
  if(! "logcounts" %in% assayNames(sce)){
    message("Logcounts assay not found; generating automatically.")
    sce <- scuttle::logNormCounts(sce)    #sce <- sclet::NormalizeData(sce)
  }
  
  # ensure the reduction method is used before setting up a cds object
  if(! reduction_method %in% names(reducedDims(sce))){
    stop("Provided SCE object lacks the specified reduction method. Please compute ", 
         reduction_method, " before running the function.")
  }
  
  # ensure the colnames of sce are not null(for running cluster_cells())
  if(is.null(colnames(sce))){
    message("The sce is missing cell names, which are required by cluster_cells; 
            generating automatically.")
    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
  }
  
  ## Initialize CellDataSet
  counts_matrix <- SummarizedExperiment::assay(sce, assay_name)
  gene_metadata <- data.frame(gene_short_name = rownames(counts_matrix), row.names = rownames(counts_matrix))
  cell_metadata <- as.data.frame(colData(sce))
  cds <- monocle3::new_cell_data_set(counts_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
  
  
  ## Add reduction result to cds object
  reducedDims(cds)[[reduction_method]] <- reducedDims(sce)[[reduction_method]]
  
  
  ## perform clustering on cells
  cds <- monocle3::cluster_cells(cds, reduction_method = reduction_method, resolution = resolution, 
                                 k = k, cluster_method = cluster_method)
  
  
  ## Trajectory and pseudotime analysis
  cds <- monocle3::learn_graph(cds)
  
  
  # a helper function to identify the root principal points
  get_earliest_principal_node <- function(cds, time_bin) {
    if (!"embryo.time.bin" %in% colnames(colData(cds))) {
      stop("The column 'embryo.time.bin' does not exist in colData(cds). Please check the input data.")
    }
    
    cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
    if (length(cell_ids) == 0) {
      stop("No cells found in the specified time_bin. Please ensure 'time_bin' is correct.")
    }
    
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
      as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
    ]
    
    return(root_pr_nodes)
  }
  
  ## pseudotime analysis(user can decide on using auto root or choosing roots on their own)
  if(use_auto_root){
    root_node <- get_earliest_principal_node(cds)
    cds <- monocle3::order_cells(cds, root_pr_nodes = root_node)
  }
  else{
    message("Please inspect the graph and choose root nodes interactively.")
    cds <- monocle3::order_cells(cds)
  }
  
  message("Monocle3 trajectory and pseudotime analysis completed.")
  return(cds)
}

