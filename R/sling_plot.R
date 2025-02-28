Genecurve_plot <- function(
    sce,
    assay_name = "logcounts",
    features,
    lineaegs = NULL,
    pseudotime.data = "slingPseudotime",
    method = c("smooth","gam"),
    point = FALSE,
    line.size = NULL,
    alpha = 0.5,
    size = 1,
    ncol = NULL,
    se = TRUE){
  
  library(ggplot2)
  library(slingshot)
  library(SingleCellExperiment)
  library(reshape2)
  
  # Ensure that pseudotime data is available
  if (!pseudotime.data %in% colnames(colData(sce))) {
    stop("The specified pseudotime data is not available in the SCE object.")
  }
  
  ## extract pseudotime data
  pseudo <- colData(sce)[[pseudotime.data]]
  pseudo <- as.data.frame(pseudo)
  
  ## extract gene expression data
  gdf <- t(assay(sce, assay_name))
  gdf <- as.data.frame(gdf[, features])
  colnames(gdf) <- features
  
  #gdf <- as.data.frame(assay(sce, assay_name)[features, ])
  #colnames(gdf) <- features
  
  df <- cbind(pseudo, gdf)
  
  ## reshape data
  d2_gene <- reshape2::melt(df, id.vars = fate_names, 
                            variable.name = "Gene", 
                            value.name = "Expression")
  
  d2 <- reshape2::melt(d2_gene, 
                       id.vars = c("Gene", "Expression"), 
                       measure.vars = fate_names, 
                       variable.name = "cell_fate", 
                       value.name = "Pseudotime")
  
  d2$Gene <- factor(d2$Gene, levels = unique(d2$Gene))
  d2 <- na.omit(d2)
  
  ## baseplot
  p <- base_plot(point, df_melt = d2, ncol)
  
  ## predict gene expression among lineages
  predictions <- data.frame()
  
  if( !is.null(lineaegs)){
    fate_names <- lineaegs
  }
  else{
    fate_names <- colnames(pseudo)
  }
  
  ## perform GSM analysis
  if (method[1] == "gam") {
    predictions <- gam_analysis_SCE(d2, features, fate_names)
  } else {
    predictions <- d2
  }
  
  p <- add_plot_layer(p, method, line.size, se = se, predictions)
  
  return(p)
}


add_plot_layer <- function(p, method, line.size, se, data_for_plot) {
  if (method[1] == "gam") {
    
    ## using prediction result to plot
    p <- p + geom_line(
      data = data_for_plot,
      aes(x = Pseudotime, y = Predicted, color = cell_fate),
      size = line.size)  
  } 
  else if (method[1] == "smooth") {
    
    # using actual gene expression data to plot
    p <- p + geom_smooth(
      data = data_for_plot,
      aes(x = Pseudotime, y = Expression, color = cell_fate),
      size = line.size, se = se)  
  }
  
  return(p)
}


base_plot <- function(
    point = FALSE, 
    df_melt, 
    alpha, 
    size, 
    ncol = NULL) {
  p <- ggplot() +
    labs(x = "Pseudotime", y = "Gene Expression", color = "Cell Fate") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(face = "bold", size = 10)) +
    facet_wrap(~ Gene, scales = "free", ncol = ncol) +
    scale_x_continuous(expand = c(0, 0))
  
  if (point) {
    p <- p + geom_point(
      data = df_melt,
      aes(Pseudotime, Expression, color = cell_fate),
      alpha = alpha, size = size,
      shape = 16)
  }
  
  return(p)
}

gam_analysis_SCE <- function(df_melt, features, fate_names) {
  
  library(mgcv)
  library(SingleCellExperiment)
  
  for (gene in features){
    for (fate in fate_names){
      # subset specific gene and lineages data from d2
      subset_df <- df_melt[df_melt$Gene == gene & df_melt$cell_fate == fate, ]
  
      # using GAM model to predict pseudotime expression
      model <- mgcv::gam(Expression ~ s(Pseudotime), data = subset_df)
  
      pseudo_seq <- seq(min(subset_df$Pseudotime), max(subset_df$Pseudotime), length.out = 200)
      grid_df <- data.frame(Pseudotime = pseudo_seq, Gene = gene, cell_fate = fate, Expression = NA)
  
      # add prediction result to final dataframe
      grid_df$Predicted <- predict(model, newdata = grid_df)
      predictions <- rbind(predictions, grid_df)
    }
  }
  predictions$Gene <- factor(predictions$Gene, levels = unique(predictions$Gene))
  
  return(predictions)
}
