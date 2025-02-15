#' Run Slingshot Trajectory Analysis
#' 
#' Runs Slingshot on a SingleCellExperiment object with specified clustering and dimensionality reduction.
#' 
#' @param sce A SingleCellExperiment object.
#' @param group The name of the column in colData(sce) containing cluster labels.
#' @param reduction_name The name of the reduced dimension method in reducedDims(sce) (e.g., "UMAP" or "PCA").
#' @param start_cluster Optional. The starting cluster for trajectory inference. Default is NULL.
#' @param end_cluster Optional. The ending cluster for trajectory inference. Default is NULL.
#' @param reverse A logical value indicating whether to reverse the pseudotime. Default is FALSE.
#' @param align_start A logical value indicating whether to align the starting pseudotime values at the maximum pseudotime. Default is FALSE.
#' @param seed Random seed for reproducibility.
#' @param show_plot Logical. If TRUE, plots the trajectories. Default is TRUE.
#' 
#' @importFrom slingshot slingPseudotime
#' @importFrom SingleCellExperiment colData reducedDims
#' @return A SingleCellExperiment object with added pseudotime information in colData.
#' @export


RunSlingshot <- function(
    sce,
    group,
    reduction_name = "UMAP",
    start_cluster = NULL,
    end_cluster = NULL,
    reverse = FALSE,
    align_start = FALSE,
    seed = 2025,
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
  
  sce <- slingshot::slingshot(data = sce, clusterLabels = group, 
                              reducedDim = reduction_name, start.clus = start_cluster, 
                              end.clus = end_cluster)
  
  
  ## reverse the pseudotime and align the start_cell to the maximum if necessary
  slingpseudotime_cols <- grep("^slingPseudotime", colnames(colData(sce)), value = TRUE)
  
  df <- as.data.frame(slingshot::slingPseudotime(sce))
  colnames(df) <- slingpseudotime_cols
  
  if (isTRUE(reverse)){
    if (isTRUE(align_start)){
      df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    }
    else{
      df <- max(df, na.rm = TRUE) - df
    }
  }
  colData(sce)[, slingpseudotime_cols] <- df
                  
  return(sce)
}

