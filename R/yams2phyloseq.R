#' Create a phyloseq Object from YAMS16 output
#'
#' @param counts unit_table.tsv count data
#' @param taxa taxa_16S_cons.tsv taxonomy data
#' @param metadata optional metadata csv/tsv
#' @param tree optional tree for Unifrac computation.
#' @export
#' @examples
#' yams2phyloseq()

yams2phyloseq <- function(counts, taxa, metadata = NULL, tree = NULL) {
  otus <- as.matrix(read.table(counts, sep="\t",row.names = 1,
                               check.names = FALSE, header = TRUE))
  taxtable <- as.matrix(read.table(taxa, sep="\t",row.names = 1,
                                   check.names = FALSE, header = TRUE))
  biom <- phyloseq(otu_table(otus,taxa_are_rows = TRUE), tax_table(taxtable))
  if (!is.null(metadata)) {
    if (!is.na(str_locate(metadata, "tsv")[1])) {
      metadata <- read.tsv(metadata)
    } else {
      metadata <- read.csv(metadata)
    }
    row.names(metadata) <- metadata[,1]
    sample_data(biom) <- metadata
  }
  if (!is.null(tree)) {
    phy_tree(biom) <- read_tree(tree)
  }
  biom
}
