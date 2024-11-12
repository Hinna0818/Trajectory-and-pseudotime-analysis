CellDimPlot <- function(
    sce,
    group,
    reduction_name = "UMAP"){
  
  ## set the cell colors
  unique_groups <- unique(colData(sce)[[group]])
  group_colors <- RColorBrewer::brewer.pal(length(unique_groups), "Set1")
  names(group_colors) <- unique_groups
  
  ## set the lineage colors
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
  
  
  ## using ggplot for visulization
  library(ggplot2)
  
  umap_1 <- colnames(reducedDims(sce)$UMAP)[1]
  umap_2 <- colnames(reducedDims(sce)$UMAP)[2]
  
  umap_data <- as.data.frame(reducedDims(sce)[[reduction_name]])
  
  set.seed(1234)
  
  p1 <- ggplot() +
    geom_point(data = umap_data, 
               aes(x = umap_1, y = umap_2, color = sce[[group]]), 
               size = 1) +  
    geom_path(data = lineage_data, aes(x = umap_1, y = umap_2, color = lineage, group = lineage),
              arrow = arrow(type = "open", length = unit(0.1, "inches")),
              size = 0.8) +
    scale_color_manual(
      values = c(group_colors, setNames(lineage_colors, paste0("Lineage ", seq_along(lineage_colors))))
    ) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(title = paste("Cells:", dim(sce)[2]))
  
  
  p2 <- p1 + guides(
    color = guide_legend(
      override.aes = list(size = 3), 
      title = "SubCellTypes and Lineages"  
    )
  ) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1) 
    ) +
    
    coord_fixed() 
  
  print(p2)
}
