#' Wrapper function for running monocle3 with default parameters.
#' 
#' @param sce SingleCellExperiment object.
#' @param assay The name of assay. Defaults to "counts".
#' @param reduction Dimensionality reduction method to use, default is "UMAP".
#' @param clusters The clustering result in SingleCellExperiment object. Defaults to NULL, in which case use Monocle clusters.
#' @param graph The name of graph slot in "metadata" of a SingleCellExperiment object. Defaults to NULL, in which case use Monocle graph is used.
#' @param partition_qval The q-value threshold for partitioning cells. Defaults to 0.05.
#' @param num_iter The number of iterations for cell clustering. Defaults to 2.
#' @param resolution The resolution parameter for cell clustering. Defaults to NULL.
#' @param cluster_method The clustering method to use. Defaults to "louvain".
#' @param k The number of nearest neighbors to consider for clustering. Defaults to 50.
#' @param use_partition Whether to use partitions to learn disjoint graph in each partition.
#' @param close_loop Whether to close loops in the graph. Defaults to TRUE.
#' @param root_pr_nodes The root nodes to order cells. User will be prompted for input before ordering cells. Defaults to NULL.
#' @param root_cells The root cells to order cells. User will be prompted for input before ordering cells. Defaults to NULL.
#' @param seed Random seed for reproducibility.
#' 
#' @importFrom igraph as_data_frame
#' @importFrom SummarizedExperiment assay assayNames
#' @export

runMonocle <- function(
    sce,
    assay = "counts",
    reduction = "UMAP",
    clusters = NULL,
    graph = NULL,
    partition_qval = 0.05,
    num_iter = 2,
    resolution = NULL,
    cluster_method = "louvain",
    k = 50,
    use_partition = TRUE,
    close_loop = TRUE,
    root_pr_nodes = NULL,
    root_cells = NULL,
    seed = 2025){

  set.seed(seed)
  
  reduction <- match.arg(reduction, c("UMAP", "tSNE", "PCA"))
  cluster_method <- match.arg(cluster_method, c("leiden", "louvain"))
  
  ## check the data before setting a cds object
  # ensure the reduction method is used before setting up a cds object
  if(! reduction %in% names(reducedDims(sce))){
    stop("The provided SCE object lacks the specified reduction method. Please compute ", 
         reduction, " before running the function.")
  }
  
  # Assign cell names if missing
  if(is.null(colnames(sce))){
    message("The sce is missing cell names, which are required by cluster_cells; 
            generating automatically.")
    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
  }
  
  ## Initialize CellDataSet (make sure the rownames of counts_matrix are gene symbols)
  counts_matrix <- SummarizedExperiment::assay(sce, assay)
  gene_metadata <- data.frame(gene_short_name = rownames(counts_matrix), row.names = rownames(counts_matrix))
  cell_metadata <- as.data.frame(colData(sce))
  cds <- monocle3::new_cell_data_set(counts_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
  
  ## check size factor
  if (!"Size_Factor" %in% colnames(cds@colData)){
    size.factor <- "size.factor"
    if (size.factor %in% colnames(colData(sce))){
      cds[["Size_Factor"]] <- colData(sce)[[size.factor]]
    }
  }
  
  ## Add reduction result and clustering result to cds object
  SingleCellExperiment::reducedDims(cds)[[reduction]] <- SingleCellExperiment::reducedDims(sce)[[reduction]]
  
  loadings <- SingleCellExperiment::reducedDims(sce)[[reduction]]
  if (length(loadings) > 0){
    slot(object = cds, name = "reduce_dim_aux")[["gene_loadings"]] <- loadings
  }
  
  stdev <- apply(reducedDims(sce)[[reduction]], 2, sd)
  if (length(stdev) > 0){
    slot(object = cds, name = "reduce_dim_aux")[["prop_var_expl"]] <- stdev
  }
  
  if (!is.null(clusters)){
    if (!is.null(graph)){
      adj_matrix <- igraph::as_adjacency_matrix(metadata(sce)[[graph]], sparse = FALSE)
      g <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_matrix, weighted = TRUE)
      
      cluster_result <- list(
        g = g,
        relations = NULL,
        distMatrix = "matrix",
        coord = NULL,
        edge_links = NULL,
        optim_res = list(
          membership = as.integer(as.factor(colData(sce)[[clusters]])),
          modularity = NA_real_
        )
      )
      if (length(unique(cluster_result$optim_res$membership)) > 1){
        cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g, cluster_result$optim_res, partition_qval)
        partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
        partitions <- as.factor(partitions)
      }
      else{
        partitions <- rep(1, ncol(sce))
      }
      names(partitions) <- colnames(cds)
      cds@clusters[["UMAP"]] <- list(
        cluster_result = cluster_result,
        partitions = partitions,
        clusters = as.factor(colData(sce)[[clusters]])
      )
      cds[["clusters"]] <- colData(sce)[[clusters]]
      cds <- monocle3:::add_citation(cds, "clusters")
      cds <- monocle3:::add_citation(cds, "partitions")
    }
    else{
      cds <- monocle3::cluster_cells(cds,
                                     reduction_method = "UMAP", partition_qval = partition_qval, k = k, cluster_method = cluster_method, 
                                     num_iter = num_iter, resolution = resolution)
      cds@clusters[["UMAP"]]$clusters <- as.factor(colData(sce)[[clusters]])
      cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters
    }
  }
  else{
    cds <- monocle3::cluster_cells(cds,
                                   reduction_method = "UMAP", partition_qval = partition_qval,
                                   k = k, cluster_method = cluster_method, num_iter = num_iter, resolution = resolution)
    cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters
  }
  
  sce[["Monocle_clusters"]] <- cds@clusters[["UMAP"]]$clusters
  sce[["Monocle_partitions"]] <- cds@clusters[["UMAP"]]$partitions
  
  
  ## learn graph
  cds <- monocle3::learn_graph(cds = cds, use_partition = use_partition, close_loop = close_loop)
  
  ## store graph results
  mst_branch_nodes <- monocle3:::branch_nodes(cds, "UMAP")
  mst_leaf_nodes <- monocle3:::leaf_nodes(cds, "UMAP")
  mst_root_nodes <- monocle3:::root_nodes(cds, "UMAP")
  pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
  
  ## pseudotime analysis
  if (is.null(root_pr_nodes) && is.null(root_cells)) {
    root_pr_nodes <- select.list(names(pps), title = "Select the root nodes to order cells, or leave blank for interactive selection:", multiple = TRUE)
    if (root_pr_nodes == "" || length(root_pr_nodes) == 0) {
      root_pr_nodes <- NULL
    }
  }
  cds <- monocle3::order_cells(cds, root_pr_nodes = root_pr_nodes, root_cells = root_cells)
  pseudotime <- cds@principal_graph_aux[["UMAP"]]$pseudotime
  pseudotime[is.infinite(pseudotime)] <- NA
  sce[["Monocle_Pseudotime"]] <- pseudotime
  
  sce@metadata$Monocle <- list(cds = cds)
  
  return(sce)
}
