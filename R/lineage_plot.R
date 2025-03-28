#' Visualize Single-Cell Trajectories on the basis of ggsc package
#' 
#' This function visualizes single-cell trajectories using UMAP (or other specified reduction method) and overlays the trajectory curves. 
#' The cells are colored by their respective lineage trajectories and the lineages are drawn as curves based on Slingshot's inferred trajectories.
#' 
#' @param sce A SingleCellExperiment object containing the results from Slingshot.
#' @param dims A vector specifying the dimensions of the reduction result to plot, default the top 2 dimensions are used.
#' @param group The name of the column in colData(sce) containing cluster labels.
#' @param reduction The name of the dimensionality reduction to plot. Default is `"UMAP"`.
#' @param lineages The lineages plotted in this function. Default all lineages are used.
#' @param cells A vector of cell names (IDs) to be included in the plot. If `NULL`, all cells are plotted.
#' @param slot The assay slot to use for dimensionality reduction, default is `"data"`. Alternatively, can be set to `"logcounts"`, `"counts"`, etc.
#' @param mapping A mapping for aesthetics such as `color`, `size`, etc., passed to `ggplot2`. By default, `NULL` and `color` will be mapped automatically.
#' @param geom A ggplot2 geometry to use for plotting. Default is `sc_geom_point`, which generates scatter points.
#' @param linewidth A numeric value indicating the line width for the trajectory paths. Default is 0.9.
#' @param ... Additional parameters passed to `ggplot2` plotting functions or the geometry.
#' 
#' @return A ggplot object displaying the UMAP plot with overlaid trajectories.
#' 
#' @examples
#' \dontrun{
#‘ pancreas <- readRDS("./pancreas_sub_sce.rds")
#’ pancreas <- RunSlingshot(sce = pancreas, group = "label", reduction = "UMAP")
#' p <- lineage_plot(sce = pancreas, group = "label", reduction = "UMAP")
#' print(p)
#' }
#' 
#' @import slingshot
#' @import SingleCellExperiment
#' @import ggplot2
#' @import ggsc
#' 
#' @export
lineage_plot <- function(
    sce,
    dims = c(1, 2),
    group,
    reduction = "UMAP",
    lineages = NULL,
    cells = NULL,
    slot = "data",
    mapping = NULL,
    geom = sc_geom_point,
    linewidth = 0.9, ...){
  
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggsc)
  
  if(is.null(lineages)){
    lineages <- names(sce@metadata$slingshot_info@curves)
  }
  
  ## set the lineage colors (adjust to ensure enough colors for the number of lineages)
  num_trajectories <- length(sce@metadata$slingshot_info@lineages)
  
  
  ## set the lineage data for the input of ggplot2
  lineage_data <- data.frame()
  
  for (i in 1:length(sce@metadata$slingshot_info@curves)) {
    temp <- as.data.frame(sce@metadata$slingshot_info@curves[[i]]$s)
    temp$lineage <- paste0("Lineage", i)
    lineage_data <- rbind(lineage_data, temp)
  }
  
  # Subset lineage_data based on the input lineage
  lineage_data <- subset(lineage_data, lineage %in% lineages)
  
  # Existing UMAP plot with labels
  p1 <- ggsc::sc_dim(object = sce, reduction = reduction, cells = cells, slot = slot,
               mapping = mapping, geom = geom)  +
    geom_path(data = lineage_data, aes(x = umap_1, y = umap_2, color = lineage, group = lineage),
              arrow = arrow(type = "open", length = unit(0.1, "inches")),
              linewidth = linewidth) 
  
  return(p1)
}
