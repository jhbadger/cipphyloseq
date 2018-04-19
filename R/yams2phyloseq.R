#' Create a phyloseq Object from YAMS16 output
#'
#' @param counts unit_table.tsv count data
#' @param taxa taxa_16S_cons.tsv taxonomy data
#' @param metadata optional metadata csv/tsv
#' @param tree optional tree for Unifrac computation.
#' @param ppm convert to parts per million
#' @export
#' @examples
#' yams2phyloseq()

yams2phyloseq <- function(counts, taxa, metadata = NULL, tree = NULL, ppm = FALSE) {
  otus <- as.matrix(read.table(counts, sep="\t",row.names = 1,
                               check.names = FALSE, header = TRUE))
  taxtable <- as.matrix(read.table(taxa, sep="\t",row.names = 1,
                                   check.names = FALSE, header = TRUE))
  biom <- phyloseq(otu_table(otus,taxa_are_rows = TRUE), tax_table(taxtable))
  if (!is.null(metadata)) {
    if (!is.na(str_locate(metadata, "csv")[1])) {
      metadata <- read.table(metadata, sep=",",row.names = 1,
                                       check.names = FALSE, header = TRUE)
    } else {
      metadata <- read.table(metadata, sep="\t",row.names = 1,
                                       check.names = FALSE, header = TRUE)
    }
    sample_data(biom) <- metadata
  }
  if (!is.null(tree)) {
    phy_tree(biom) <- read_tree(tree)
  }
  if (ppm) {
    biom <- transform_sample_counts(biom, function(x)
      round(1E6 * x/sum(x)))
  }
  biom
}

#' Load QIIME2 data into metagenomeSeq
#'
#' @param tableqza table.qza qiime2 count data
#' @param taxonomyqza taxonomy.qza qiime2 taxonomy data
#' @param metadata metadata csv/tsv
#' @export
#' @examples
#' qiime2MR()

qiime2MR <- function(tableqza, taxonomyqza, metadata, ppm=TRUE) {
  tmpd <- tempfile()
  if (!dir.exists(tmpd)) {
    dir.create(tmpd)
  }
  system(str_c("unzip -d ", tmpd, " -j ",
    tableqza, " *.biom"),
    ignore.stdout = TRUE, ignore.stderr = TRUE)
  biom <- read_biom(str_c(tmpd, "/feature-table.biom"))
  biom <- biom2MRexperiment(biom)
  if (ppm) {
    counts <- MRcounts(biom)
    counts <- round(sweep(counts, 2, colSums(counts), FUN="/")*1e6)
    biom <- newMRexperiment(counts=counts)
  }
  system(str_c("unzip -d ", tmpd, " -j ",taxonomyqza, " *.tsv"),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  taxonomy <- read.table(str_c(tmpd, "/taxonomy.tsv"), sep="\t", header=TRUE, row.names = 1)
  system(str_c("rm -rf ", tmpd))
  ts <- data.frame(str_split_fixed(taxonomy$Taxon, ';', n=8)[,1:7])
  row.names(ts) <- rownames(taxonomy)
  colnames(ts) <- c("Rank0","Rank1","Rank2","Rank3","Rank4","Rank5","Rank6")
  ts$LKT <- apply(ts, 1, function(x) {
    y=str_replace_all(x,str_c("uncultured bacterium|uncultured organism|", "
                              uncultured$|unidentified|Incertae Sedis"),"")
    ifelse(length(which(str_length(y) < 6))>0,
           y[min(which(str_length(y) < 6)) - 1], y[7])})

  fData(biom) <- ts
  metadata <- read.table(metadata,
                         sep = ifelse(is.na(
                           str_locate(metadata, "csv")[1]),
                           "\t",","),
                         row.names = 1,
                         check.names = FALSE,
                         header = TRUE)
  pData(biom) <- metadata
  biom
}

#' Convert metagenomeSeq to phyloseq
#'
#' @param mrexp metagenomeSeq experiment object
#' @param unrootedtreeqza optional unrootedtree.qza for Unifrac
#' @export
#' @examples
#' mr2phyloseq()

mr2phyloseq <- function(mrexp, unrootedtreeqza=NULL) {
  counts <- MRcounts(mrexp)
  metadata <- pData(mrexp)
  taxonomy <- fData(mrexp)
  phy <- phyloseq(otu_table(counts, taxa_are_rows = TRUE),
                  sample_data(metadata))
  tax_table(phy) <- as.matrix(taxonomy[taxa_names(phy),])
  if (!is.null(unrootedtreeqza)) {
    tmpd <- tempfile()
    if (!dir.exists(tmpd)) {
      dir.create(tmpd)
    }
    system(str_c("unzip -d ", tmpd, " -j ",
    unrootedtreeqza, " *.nwk"),
    ignore.stdout = TRUE, ignore.stderr = TRUE)
    phy_tree(phy) <- read_tree(str_c(tmpd,"/tree.nwk"))
    system(str_c("rm -rf ", tmpd))
  }
  phy
}
