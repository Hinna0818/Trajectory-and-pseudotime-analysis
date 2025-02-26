#' Visualize Single-Cell Trajectories on the basis of ggsc package
#' 
#' This function visualizes single-cell trajectories using UMAP (or other specified reduction method) and overlays the trajectory curves. 
#' The cells are colored by their respective lineage trajectories and the lineages are drawn as curves based on Slingshot's inferred trajectories.
#' 
#' @param sce A SingleCellExperiment object containing the results from Slingshot.
#' @param dims A vector specifying the dimensions of the reduction result to plot, default the top 2 dimensions are used.
#' @param The name of the column in colData(sce) containing cluster labels.
#' @param reduction The name of the dimensionality reduction to plot. Default is `"UMAP"`, but it can be changed to other reductions like `"PCA"` or `"tSNE"`.
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
#' # Assuming `sce` is a SingleCellExperiment object with Slingshot results
#' p <- CDP(sce = sce, group = "label", reduction = "UMAP")
#' print(p)
#' }
#' 
#' @seealso \code{\link{sc_dim}}, \code{\link{sc_geom_point}}, \code{\link{slingshot}}
#' 
#' @import slingshot
#' @import SingleCellExperiment
#' @import ggplot2
#' @import ggsc
#' 
#' @export
CDP <- function(
    sce,
    dims = c(1, 2),
    group,
    reduction = "UMAP",
    cells = NULL,
    slot = "data",
    mapping = NULL,
    geom = sc_geom_point,
    linewidth = 0.9, ...){
  
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggsc)
  
  ## set the lineage colors (adjust to ensure enough colors for the number of lineages)
  num_trajectories <- length(sce@metadata$slingshot_info@lineages)
  
  lineage_colors <- RColorBrewer::brewer.pal(min(num_trajectories, 8), "Set1")
  if (num_trajectories > 8) {
    lineage_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(num_trajectories)
  }
  
  ## set the lineage data for the input of ggplot2
  lineage_data <- data.frame()
  
  for (i in 1:length(sce@metadata$slingshot_info@curves)) {
    temp <- as.data.frame(sce@metadata$slingshot_info@curves[[i]]$s)
    temp$lineage <- paste0("Lineage ", i)
    lineage_data <- rbind(lineage_data, temp)
  }
  
  # Existing UMAP plot with labels
  p1 <- sc_dim(object = sce, reduction = reduction, cells = cells, slot = slot,
               mapping = mapping, geom = geom)  +
    geom_path(data = lineage_data, aes(x = umap_1, y = umap_2, color = lineage, group = lineage),
              arrow = arrow(type = "open", length = unit(0.1, "inches")),
              linewidth = linewidth) 
  
  return(p1)
}
