#' Make a heat map of a biom object with clustering
#'
#' @param biom phyloseq object
#' @param rank string giving taxonomic rank to compare at
#' @param group vector identifying the group
#' @param title string giving the title
#' @param ntop maximum number of rows
#' @export
#' @examples
#' heatmap_biom()

heatmap_biom <- function(biom, rank="LKT", group,
                          title="", ntop=40) {
  biom <- tax_glom(biom, rank)
  taxa_names(biom) <- make.names(tax_table(biom)[,rank],
                                     unique = TRUE)
  ntop <- min(ntaxa(biom)/2, ntop)
  obj <- phyloseq_to_metagenomeSeq(biom)
  trials <- pData(obj)[ , which(colnames(pData(obj))==group)]
  heatmapColColors <- brewer.pal(12,"Set3")[as.integer(factor(trials))]
  heatmapCols <- colorRampPalette((brewer.pal(9, "YlOrRd")))(50)
  par(oma=c(5,4,3,12)+0.1, xpd=TRUE)
  plotMRheatmap(obj, n = ntop, norm = TRUE, log = TRUE,
                fun = sd, cexRow = 0.5, cexCol = 0.5, trace = "none",
                col = heatmapCols, ColSideColors = heatmapColColors,
                key=FALSE, main=title, cex.main=0.35)
  cores <- unique(heatmapColColors)
  labels <- unique(trials)
  legend("topright", inset=c(1,0), legend=labels, fill=cores,
         border=TRUE, bty="n", y.intersp = 0.7, cex=0.5)
}
