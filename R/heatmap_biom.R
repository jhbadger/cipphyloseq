#' Make a heat map of a biom object with clustering
#'
#' @param biom phyloseq object
#' @param rank string giving taxonomic rank to compare at
#' @param group  vector identifying the group for each slice
#' @param glom run tax_glom on biom at rank
#' @param norm normalize (via metagenomeSeq) data
#' @param log log-transform data
#' @param n maximum number of rows to include
#' @cexRow size of row font
#' @cexCol size of column font
#' @reorder reorder columns by group
#' @minVar minimum variance to consider row
#' @export
#' @examples
#' heatmap_biom()

heatmap_biom <- function(biom, rank, group, glom=TRUE, reorder = TRUE,
                         title = "", lowCol = "blue", highCol="tomato",
                         log = TRUE) {
  if (glom) {
    biom <- tax_glom(biom, rank)
    taxa_names(biom) <- make.unique(as.data.frame(tax_table(biom))[[rank]] %>% as.character)
  }
  counts <- as.matrix(otu_table(biom))
  colnames <- colnames(counts)
  cdf <- data.frame(scale(counts))
  if (log) {
    cdf <- log(cdf + 1)
  }
  colnames(cdf) <- colnames
  cdf$vars <- rowVars(as.matrix(cdf))
  cdf$taxa <- rownames(cdf)
  cdf_melt <- melt(cdf, id=c("taxa", "vars"), variable_name = "sample")
  cdf_melt <- cdf_melt[order(cdf_melt$vars, decreasing = TRUE),]
  if (reorder) {
    
  }
  ggplot(data = cdf_melt, aes(x = sample, y = taxa)) + 
    geom_tile(aes(fill = value), size = 1) + 
    scale_fill_gradient(low = lowCol, high = highCol) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, ),) 
    ggtitle(title)
}
