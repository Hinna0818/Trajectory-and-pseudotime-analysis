#' Generate gene expression curves across pseudotime
#'
#' This function visualizes the expression of specified genes across pseudotime using either 
#' smoothed line plots or GAM (Generalized Additive Models) predictions. It generates a 
#' ggplot showing the expression of selected genes along pseudotime for different cell fates.
#'
#' @param sce A SingleCellExperiment object containing the expression data and pseudotime.
#' @param assay_name The assay to use for extracting gene expression data (default: "logcounts").
#' @param features A vector of gene names for which to plot expression curves.
#' @param lineaegs A vector of cell fates to visualize (optional). If NULL, all available fates will be used.
#' @param pseudotime.data The column name in the SingleCellExperiment object containing pseudotime data (default: "slingPseudotime").
#' @param method A character vector specifying the method for curve fitting. Options are "smooth" (default) or "gam" (Generalized Additive Model).
#' @param point A logical value indicating whether to include points for individual cells (default: FALSE).
#' @param line.size Line thickness for the plot (default: NULL).
#' @param alpha Transparency for the points in the plot (default: 0.5).
#' @param size Size of the points in the plot (default: 1).
#' @param ncol Number of columns to use in facet_wrap (default: NULL).
#' @param se A logical value indicating whether to include the standard error for smoothing (default: TRUE).
#'
#' @return A ggplot object displaying gene expression curves across pseudotime.
#' @export
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
  
  if(! "logcounts" %in% assayNames(sce)){
    sce <- NormalizeData(sce)
  }
  
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



#' Add plot layers based on method choice
#'
#' This function adds the appropriate plot layer to the base plot depending on the selected method. 
#' It either adds a GAM-based line plot or a smoothed curve for the gene expression across pseudotime.
#'
#' @param p The base ggplot object.
#' @param method A character vector specifying the plot method, either "smooth" or "gam".
#' @param line.size Size of the plot line (default: NULL).
#' @param se A logical value indicating whether to display the standard error (default: TRUE).
#' @param data_for_plot The data used for plotting, either predictions or raw data.
#'
#' @return A ggplot object with the added plot layer.
#' @export
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


#' Create the base plot for gene expression over pseudotime
#'
#' This function generates the base plot framework, which includes labels, faceting by gene, and 
#' basic styling. It can optionally add points to represent individual cell data.
#'
#' @param point A logical value indicating whether to add points for individual cells (default: FALSE).
#' @param df_melt A data frame containing the melted version of the gene expression data.
#' @param alpha The transparency level for points (default: 0.5).
#' @param size The size of the points (default: 1).
#' @param ncol The number of columns for faceting (default: NULL).
#'
#' @return A ggplot object representing the base plot.
#' @export
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


#' Perform GAM analysis on gene expression data
#'
#' This function fits a Generalized Additive Model (GAM) to the gene expression data for each
#' gene and lineage, predicting gene expression across pseudotime.
#'
#' @param df_melt A data frame containing melted gene expression data, including pseudotime and cell fate.
#' @param features A vector of gene names to analyze.
#' @param fate_names A vector of cell fate names to consider for each gene.
#'
#' @return A data frame with predicted gene expression values based on the fitted GAM models.
#' @export
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
