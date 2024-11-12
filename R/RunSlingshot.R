#' Run Slingshot Trajectory Analysis
#' 
#' Runs Slingshot on a SingleCellExperiment object with specified clustering and dimensionality reduction.
#' 
#' @param sce A SingleCellExperiment object.
#' @param group The name of the column in colData(sce) containing cluster labels.
#' @param reduction_name The name of the reduced dimension method in reducedDims(sce) (e.g., "UMAP" or "PCA").
#' @param start_cluster Optional. The starting cluster for trajectory inference.
#' @param end_cluster Optional. The ending cluster for trajectory inference.
#' @param dims A vector of integers specifying which dimensions to use. Default is the first two dimensions.
#' @param seed Random seed for reproducibility.
#' @param show_plot Logical. If TRUE, plots the trajectories. Default is TRUE.
#' 
#' @importFrom slingshot slingshot slingPseudotime slingBranchID
#' @importFrom SingleCellExperiment colData reducedDims
#' @importFrom graphics plot lines
#' @return A SingleCellExperiment object with added pseudotime information in colData.
#' @export


RunSlingshot <- function(
    sce,
    group,
    reduction_name = "UMAP",
    start_cluster = NULL,
    end_cluster = NULL,
    dims = 1:2,
    seed = 2024-11-08,
    show_plot = TRUE){
  
  reductition_name <- match.arg(reduction_name, c("UMAP", "PCA", "tSNE"))
  
  ## check the type of data
  if (!is(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  ## check the cluster result in sce obj
  if (!group %in% colnames(colData(sce))) {
    stop("The specified cluster column does not exist in colData(sce).")
  }
  
  ## check the reduction result in sce obj
  if (!reduction_name %in% names(reducedDims(sce))) {
    stop("The specified reduction method does not exist in reducedDims(sce).")
  }
  
  
  ## run slingshot
  set.seed(seed)
  
  sds <- slingshot::slingshot(data = sce, clusterLabels = group, 
                              reducedDim = reduction_name, start.clus = start_cluster, 
                              end.clus = end_cluster)
  
  ssds <- slingshot::SlingshotDataSet(sds) 
  
  # Store Slingshot results
  colData(sce)$slingPseudotime <- as.data.frame(slingshot::slingPseudotime(ssds))
  colData(sce)$slingBranchID <- slingshot::slingBranchID(ssds)
  metadata(sce)$slingshot_info <- ssds
  
  return(sce)
}
  
