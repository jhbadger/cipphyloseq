#' Make a scatter PCA/PCOA plot of phyloseq object
#'
#' @param biom phyloseq object
#' @param method data reduction method to use (default PCoA; see phyloseq ordinate)
#' @param distance distance measure to use (default unifrac)
#' @param color_category variable in sample data (if any) to color points by
#' @param shape_category variable in sample data (if any) to distinguish points by symbol
#' @param weighted use weighted distance (normally only useful for unifrac)
#' @param title title of plot
#' @param color vector of colors to use if you don't like the defaults
#' @param shapes vector of shapes to use if you don't like the defaults
#' @export
#' @examples
#' pcoa_biom()

pcoa_biom <- function(biom, method = "PCoA", distance = "unifrac", color_category = NULL,
                      shape_category=NULL, weighted = FALSE, title = NULL, colors = NULL,
                      shapes = NULL) {
  ordu <- ordinate(biom, method = method, distance = distance, weighted = weighted)
  p <- plot_ordination(biom, ordu, color = color_category, shape = shape_category)
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values=colors)
  }
  if(!is.null(shapes)) {
    p <- p + scale_shape_manual(values=shapes)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p +  theme_light()
}

#' Make a scatter PCA/PCOA plot of phyloseq object with complete biom and subset of samples -- useful for conserved axes
#' 
#' @param biom phyloseq object
#' @param method data reduction method to use (default PCoA; see phyloseq ordinate)
#' @param distance distance measure to use (default unifrac)
#' @param color_category variable in sample data (if any) to color points by
#' @param shape_category variable in sample data (if any) to distinguish points by symbol
#' @param weighted use weighted distance (normally only useful for unifrac)
#' @param title title of plot
#' @param color vector of colors to use if you don't like the defaults
#' @param shapes vector of shapes to use if you don't like the defaults
#' @export
#' @examples
#' pcoa_subset_biom

pcoa_subset_biom <- function(biom, method = "PCoA", distance = "unifrac", samples, 
                             color_category = NULL, shape_category=NULL, weighted = FALSE, 
                             title = NULL, colors = NULL, shapes = NULL, seed = 100) {
  set.seed(seed)
  ordu <- ordinate(biom, method = method, distance = distance, weighted = weighted)
  plt <- plot_ordination(biom, ordu, color = color_category, 
                         shape = shape_category, justDF = TRUE)
  minX <- min(plt$Axis.1)*1.1
  minY <- min(plt$Axis.2)*1.1
  maxX <- max(plt$Axis.1)*1.1
  maxY <- max(plt$Axis.2)*1.1
  plt <- plt[samples, ]
  plt <- plt[complete.cases(plt),]
  p <- ggplot(plt, aes_string(x="Axis.1", y="Axis.2", color=color_category,
                              shape = shape_category))+geom_point()
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values=colors)
  }
  if(!is.null(shapes)) {
    p <- p + scale_shape_manual(values=shapes)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + theme_light()
  p <- p + coord_cartesian(xlim = c(minX, maxX),
                           ylim = c(minY, maxY)) 
  p <- p + xlab(paste0(paste0("Axis.1 [",
                       round(ordu$values$Relative_eig[1]*100,1)), "%]"))
  p + 
    ylab(paste0(paste0("Axis.2 [",
                       round(ordu$values$Relative_eig[2]*100,1)), "%]"))
}
