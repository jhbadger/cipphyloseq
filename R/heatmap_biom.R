#' Make a heat map of a biom object with clustering
#'
#' @param biom phyloseq object
#' @param rank string giving taxonomic rank to compare at
#' @param group  vector identifying the group for each slice
#' @export
#' @examples
#' heatmap_biom()

heatmap_biom <- function(biom, rank, group, title = "") {
  biom <- tax_glom(biom, rank)
  plot_heatmap(biom, taxa.label = rank, sample.label = group) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     size = min(8, 350/nsamples(biom))), 
          axis.text.y = element_text(size = min(10, 350/ntaxa(biom))))  +
    ggtitle(title) 
}
