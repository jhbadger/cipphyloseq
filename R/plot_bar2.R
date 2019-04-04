
#' More flexible version of phyloseq's plot_bar
#'
#'
#' @param biom a phyloseq object.
#' @param fill string giving taxonomic level (from tax_table) to plot at.
#' @param title title of plot
#' @param colors vector of color names to use
#' @param minPercent optional limit to exclude low abundant taxa
#' @param pickTaxa optional list of taxa to limit graph to.
#' @param excludeTaxa optional list of taxa to exclude
#' @param border_color border color of bars
#' @param log log transform data
#' @param ppm translate data into parts per million rather than percent
#' @param rbiom return a phyloseq data object
#' @param showPercents limit graph to c(min,max) percent
#' @export
#' @examples
#' plot_bar2()

plot_bar2 <- function (biom, fill = NULL, title = NULL, colors = NULL, minPercent = NULL,
                             pickTaxa = NULL, excludeTaxa = NULL, border_color = "black",
                             log = FALSE, ppm = FALSE, rbiom = FALSE, showPercents = NULL,
                       nomisc=FALSE)
{
  units <- 100
  ylabel <- "Relative Abundance (%)"
  if (ppm) {
    units <- 1e+06
    ylabel <- "Relative Abundance (ppm)"
    if (!is.null(minPercent)) {
      minPercent <- minPercent * 10000
    }
  }
  physeq <- transform_sample_counts(biom, function(x) units *
                                      x/sum(x))
  if (!is.null(pickTaxa)) {
    otus <- tax_table(physeq) %>% as.data.frame
    otus <- row.names(otus)[otus[, fill] %in% pickTaxa]
    physeq <- prune_taxa(otus, physeq)
    physeq <- transform_sample_counts(physeq, function(x) units *
                                        x/sum(x))
  }
  if (!is.null(excludeTaxa)) {
    otus <- tax_table(physeq) %>% as.data.frame
    otus <- row.names(otus)[!(otus[, fill] %in% excludeTaxa)]
    physeq <- prune_taxa(otus, physeq)
    physeq <- transform_sample_counts(physeq, function(x) units *
                                        x/sum(x))
  }
  if (log) {
    otu_table(physeq)[otu_table(physeq) == 0] = 1
    otu_table(physeq) <- log(otu_table(physeq))
    physeq <- transform_sample_counts(physeq, function(x) log(units) *
                                        x/sum(x))
  }
  physeq <- tax_glom(physeq, fill)
  if (!is.null(minPercent)) {
    physeq <- filter_taxa(physeq, function(x) mean(x) > minPercent,
                         TRUE)
    physeq <- transform_sample_counts(physeq, function(x) round(x))
    misc_counts <- units-colSums(otu_table(physeq))
    new_table <- rbind(otu_table(physeq), Misc_Low_Abundance=misc_counts)
    if (!nomisc) {
      new_tax <- rbind(as.data.frame(tax_table(physeq), stringsAsFactors = FALSE),
                       Misc_Low_Abundance = rep("Misc_Low_Abundance", ncol(tax_table(physeq))))
      physeq <- phyloseq(otu_table(new_table, taxa_are_rows = TRUE),
                         tax_table(as.matrix(new_tax)))
    } else {
      physeq <- transform_sample_counts(physeq, function(x) units *
                                          x/sum(x))
    }
  }
  if (rbiom) {
    return(physeq)
  }
  mdf <- psmelt(physeq)
  fill_sum <- aggregate(as.formula(paste("Abundance ~", fill)),
                        mdf, sum)
  mdf[, fill] <- factor(mdf[, fill], levels = fill_sum[order(-fill_sum$Abundance),
                                                       fill])
  mdf$Sample <- factor(mdf$Sample, levels = sample_names(physeq))
  p <- ggplot(mdf, aes_string(x = "Sample", y = "Abundance",
                              fill = fill))
  p <- p + geom_bar(stat = "identity", position = "stack",
                    color = border_color, size = 0)
  p <- p + ylab(ylabel)
  p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                 axis.ticks.x = element_blank())
  p <- p + theme(plot.margin = unit(c(0, 1, 1, 0.5), "lines"))
  p <- p + theme(strip.background = element_rect(fill = "white",
                                                 colour = NA))
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  if (!is_null(showPercents)) {
    p <- p + coord_cartesian(ylim = c(showPercents[1], showPercents[2]))
  }
  p <- p + theme(panel.background = element_rect(fill = "white",
                                                 colour = "white"))
  p <- p + theme(plot.margin = margin(c(1, 1, 3, 1), unit = "lines"))
  p <- p + theme(axis.line = element_blank())
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
  if (!is.infinite(x)) {
    xend <- max(which(categories == matching))
    my_bracket(plot, text, x, xend, y, cex)
  } else plot
}
