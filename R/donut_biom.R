#' Make a donut graph of two taxonomic levels
#'
#' @param biom phyloseq object
#' @param group  vector identifying the group for each slice
#' @param labels vector of labels for individual slices
#' @param inner_colors colors four inner circle
#' @param outer_colors colors for outer circle
#' @param threshold threshold percent below which don't display
#' @param title title of graph
#' @export
#' @examples
#' donut_biom()

donut_biom <- function(biom, samples=NULL, rank1, rank2, inner_colors,
                       outer_colors, threshold = 1, title="") {

  if (!is.null(samples)) {
    biom <- prune_samples(samples, biom)
  }
  biom <- tax_glom(biom, rank2)

  counts <- data.frame(weight = as.numeric(otu_table(biom)),
                       rank1 = as.data.frame(tax_table(biom))[[rank1]],
                       Taxon = as.data.frame(tax_table(biom))[[rank2]])
  counts$weight <- 100.0*counts$weight/sum(counts$weight)
  counts <- filter(counts, weight>=threshold)
  counts$weight <- 100.0*counts$weight/sum(counts$weight)
  counts <- counts %>% group_by(rank1,Taxon) %>% summarise(weight = sum(weight))
  counts$fraction <- counts$weight / sum(counts$weight)
  counts <- counts[order(counts$fraction), ]
  counts$ymax <- cumsum(counts$fraction)
  counts$ymin <- c(0, head(counts$ymax, n=-1))
  counts$yavg <- (counts$ymax+counts$ymin)/2
  colors <- append(inner_colors[1:length(unique(counts$rank1))],
                   outer_colors[1:length(unique(counts$Taxon))])
  p <- ggplot(counts) + 
    geom_rect(aes(ymin=ymin, ymax=ymax, xmax=4, xmin=3, fill=Taxon, color=Taxon)) + 
    scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
    geom_rect(aes(ymin=ymin, ymax=ymax, xmax=3, xmin=0, fill=rank1, color=rank1)) +
    xlim(c(0, 4)) + 
    theme(aspect.ratio=1) +
    coord_polar(theta="y") 
  p <- p + theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(axis.line = element_line(colour = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  p <- p + ggtitle(title) + xlab("") + ylab("")
  p
}


