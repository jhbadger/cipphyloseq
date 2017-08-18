#' Make a donut graph of two taxonomic levels
#'
#' @param biom a phyloseq object.
#' @param taxlevel1 string giving taxonomic level (from tax_table) to plot inner pie at.
#' @param taxlevel2 string giving taxonomic level (from tax_table) to plot outer pie at.
#' @param variable variable in sample data to restrict samples to.
#' @param condition value of variable to restrict samples to
#' @param threshold threshold percentage above which will be included in graph
#' @export
#' @examples
#' donut_biom
#'
donut_biom <- function(biom, taxlevel1, taxlevel2, variable, condition, threshold = 0) {
  library(ggiraphExtra)
  library(ggiraph)
  library(phyloseq)
  t <- tax_glom(biom, taxlevel2)
  counts <- as.data.frame(tax_table(t))[,c(taxlevel1, taxlevel2)]
  t <- prune_samples(sample_names(t)[sample_data(t)[[variable]]==condition], t)
  counts$count <- as.numeric(otu_table(merge_samples(t, variable)))
  counts$count <- 100.0*counts$count / sum(counts$count)
  counts <- counts[counts$count>threshold,]
  colnames(counts) <- c("R1","R2","count")

  ggPieDonut(counts, aes(pies=R1, donuts = R2, count = count))
}
