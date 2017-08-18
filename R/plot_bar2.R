
#' More flexible version of phyloseq's plot_bar
#'
#'
#' @param biom a phyloseq object.
#' @param fill string giving taxonomic level (from tax_table) to plot at.
#' @param title title of plot
#' @param colors vector of color names to use
#' @param maxTaxa optional limit to number of taxa to display
#' @export
#' @examples
#' plot_bar2()

plot_bar2 <- function (biom, fill = NULL, title = NULL, colors = NULL, maxTaxa = NULL)
{
  library(phyloseq)
  library(ggplot2)
  physeq <- transform_sample_counts(biom, function(x) 100*x / sum(x))
  mdf <- psmelt(tax_glom(physeq, fill))
  fill_sum <- aggregate(as.formula(paste("Abundance ~",fill)), mdf, sum)
  mdf[,fill] <- factor(mdf[,fill], levels=fill_sum[order(-fill_sum$Abundance), fill])
  if (!is.null(maxTaxa)) {
    physeq <- prune_taxa(levels(mdf[,fill])[1:maxTaxa], physeq)
    physeq <- transform_sample_counts(physeq, function(x) 100*x / sum(x))
    mdf <- psmelt(tax_glom(physeq, fill))
    fill_sum <- aggregate(as.formula(paste("Abundance ~",fill)), mdf, sum)
    mdf[,fill] <- factor(mdf[,fill], levels=fill_sum[order(-fill_sum$Abundance), fill])
  }
  mdf$Sample <- factor(mdf$Sample, levels=sample_names(physeq))
  p <- ggplot(mdf, aes_string(x = "Sample", y = "Abundance", fill = fill))
  p <- p + geom_bar(stat = "identity", position = "stack", color = "black", size = 0)
  p <- p + ylab("Relative Abundance (%)")
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),  axis.ticks.x=element_blank())
  p <- p + theme(plot.margin = unit(c(0, 1, 1, 0.5), "lines"))
  p <- p + theme(strip.background = element_rect(fill = "white", colour = NA))
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values=colors)
  }
  p <- p + theme(panel.background = element_rect(fill = 'white', colour = 'white'))
  return(p)
}

#' Adds a bracket with a label to annotate a ggplot2 plot
#'
#' For example, you can label groups of samples in a large bar chart.
#' More typically you will use auto_bracket() which automatically sets
#' the x and xend based on a category label.
#'
#' @param plot a ggplot2 plot.
#' @param text the text to label the bracket with.
#' @param x the x coordinate of the start of the bracket.
#' @param xend the x coordinate of the end of the bracket.
#' @param y the y coordinate of the bracket.
#' @param cex the size of the text (default 2.2)
#' @export
#' @examples
#' my_bracket()

my_bracket <- function(plot, text, x, xend, y, cex = 2.2) {
  plot + annotate("segment", x = x, xend = xend, y = y, yend = y) +
    annotate("text", label = text, x = x+(xend-x)/2, y= y-3, cex=cex)
}

#' Adds a bracket automatically with a label to annotate a ggplot2 plot
#'
#' For example, you can label groups of samples in a large bar chart.
#' auto_bracket automatically positions the bracket assuming your list
#' of categories is in the same order as your samples (which they would
#' if they were part of a metadata sheet, or a sample_data in phyloseq)
#'
#' @param plot a ggplot2 plot.
#' @param matching the category value for the bracket
#' @param categories the vector of categories (one per sample)
#' @param y the y coordinate of the bracket.
#' @param cex the size of the text (default 2.2)
#' @param text the value of the text (if it isn't the same as matching)
#' @export
#' @examples
#' auto_bracket()

auto_bracket <- function(plot, matching, y, categories, cex = 2.2, text = matching) {
  x <- min(which(categories == matching))
  xend <- max(which(categories == matching))
  my_bracket(plot, text, x, xend, y, cex)
}
