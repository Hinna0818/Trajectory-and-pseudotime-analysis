#' Visualize pseudotime data on reduction plots
#'
#' @param sce A SingleCellExperiment object after runSlingshot which contains the pseudotime and reduction data.
#' @param reduction A string specifying the reduction method to be used for plotting. Default is `NULL`.
#' @param pseudotime The name of the column in the colData of `sce` that contains the pseudotime data. Default is `"slingPseudotime"`.
#' @param ... Additional parameters to be passed to `ggplot2` and other plot settings.
#'
#' @return A `ggplot` object representing the pseudotime visualization.
#'
#' @examples
#' p <- pseudo_plot(sce = sce, reduction = "UMAP")
#' print(p)
#'
#' @import ggplot2
#' @import SingleCellExperiment
#' @importFrom tidyr pivot_longer
#' @importFrom viridis scale_color_viridis_c
#'
#' @export
pseudo_plot <- function(
    sce,
    reduction = NULL,
    pseudotime = "slingPseudotime",
    ...){
  
  library(SingleCellExperiment)
  library(ggplot2)
  
  reduction <- match.arg(reduction, c("UMAP", "PCA", "tSNE"))
  
  if( !is(sce, "SingleCellExperiment")){
    stop("The input must be a SingleCellExperiment object.")
  }
  
  if( !pseudotime %in% colnames(colData(sce))){
    stop("SCE lacks pseudotime data. Please run 'runSlingshot' before this function.")
  }
  
  
  ## extract pseudotime data
  sling <- sce$slingPseudotime
  sling$cell <- rownames(sling)
  d2 <- tidyr::pivot_longer(sling, 
                            cols = starts_with("Lineage"),
                            names_to = "Lineage", 
                            values_to = "Pseudotime")
  
  ## extract reduction data
  reduction_data <- as.data.frame(reducedDims(sce)[[reduction]])
  reduction_data$cell <- rownames(reduction_data)
  
  ## merge pseudotime data and reduction data
  merged <- merge(d2, reduction_data, by = "cell")
  
  ## use ggplot2 to plot psedotime data of cells
  p <- ggplot(merged, aes(x = umap_1, y = umap_2, color = Pseudotime)) +
    geom_point() +
    facet_wrap(~Lineage, ncol = 3) +  
    scale_color_viridis_c(option = "C") + 
    theme_bw() +
    theme(
      axis.text = element_blank(),    
      axis.title = element_blank(),   
      axis.ticks = element_blank(),   
      panel.grid = element_blank(),   
      strip.text = element_text(size = 12, face = "bold"),  
      legend.position = "right",      
      legend.title = element_text(size = 12),  
      legend.text = element_text(size = 10)    
    ) +
    labs(
      title = NULL,  
      x = NULL,      
      y = NULL       
    )
  
  return(p)
}