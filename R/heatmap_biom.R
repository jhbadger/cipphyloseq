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

heatmap_biom <- function(biom, rank, group, glom=TRUE, norm = TRUE, log=TRUE, n = ntaxa(biom), cexRow=0.5, cexCol=0.5, reorder=TRUE, minVar = 2) {
  if (glom) {
    biom <- tax_glom(biom, rank)
    taxa_names(biom) <- make.unique(as.data.frame(tax_table(biom))[[rank]] %>% as.character)
  }
  goodtaxa <- taxa_names(biom)[apply(as.data.frame(otu_table(biom)), 1, var) >= minVar]
  biom <- prune_taxa(goodtaxa, biom)
  group_num <- as.character(as.numeric(sample_data(biom)[[group]]))
  plotMRheatmap(obj = phyloseq_to_metagenomeSeq(biom), norm=norm, log=log, n = n, trace = "none",
                cexRow=cexRow, cexCol=cexCol, key=TRUE, keysize=0.5, Colv = reorder,
                ColSideColors=group_num, key.title=NA, key.ylab=NA, key.xlab=NA,
                margins=c(8,8), fun=mad)
  par(lend = 1, xpd=TRUE)
  legend(x=0,y=2.1, legend=unique(sample_data(biom)[[group]]),
         col=unique(sample_data(biom)[[group]]), lty=1, lwd=10)
}
