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
