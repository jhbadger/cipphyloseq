#' Make a donut graph of two taxonomic levels
#'
#' biom
#' group  vector identifying the group for each slice
#' labels vector of labels for individual slices
#' col    colors for each group
#' radius radius for inner and outer pie (usually in [0,1])
#' threshold threshold percent below which don't display
#' main title of graph
#' donut_biom

donut_biom <- function(biom, samples, rank1, rank2, inner_colors,
                       outer_colors, radius = c(.7, 1), threshold = 1, main="") {

  if (!is.null(samples)) {
    biom <- prune_samples(samples, biom)
  }
  biom <- tax_glom(biom, rank2)
  biom <- merge_samples(biom, rep(TRUE, length(sample_names(biom))), sum)
  biom <- transform_sample_counts(biom, function(x) 100*x / sum(x))
  biom <- prune_taxa(taxa_names(biom)[otu_table(biom)>=threshold], biom)
  outer_labels <- as.data.frame(tax_table(biom))[[rank2]]
  x <- as.numeric(otu_table(biom))
  plot.new()
  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L],
      col = outer_colors, labels = outer_labels, main = main)
  biom <- tax_glom(biom, rank1)
  x <- as.numeric(otu_table(biom))
  inner_labels <- as.data.frame(tax_table(biom))[[rank1]]
  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L],
      col = inner_colors, labels=NA)
  legend("topright", legend=inner_labels, cex=0.8, fill=inner_colors)
  par(mar=c(5,3,2,2)+0.1)
}

