#' Run DESeq2 on a phyloseq object and return tibble of significantly abundant OTUs
#'
#' @param biom a phyloseq object.
#' @param taxlevel string giving taxonomic level (from tax_table) to compare at.
#' @param group string giving column in sample data to compare
#' @param alpha p-value cutoff
#' @param threshold minimum percentage of OTU to report
#' @param glom run tax_glom on biom to get to desired level (FALSE if already at right level)
#' @export
#' @examples
#' deseq2_biom()

deseq2_biom <- function(biom, taxlevel, group, alpha = 0.01, threshold=1, glom=TRUE,
                        pseudocounts = 0, maxgroups = 5) {
  library(phyloseq)
  library(DESeq2)
  if (glom) {
    biom <- tax_glom(biom, taxrank = taxlevel)
  }
  norm <-  transform_sample_counts(biom, function(x) 100*x / sum(x))
  biom <- prune_taxa(taxa_names(norm)[rowMeans(as.data.frame(otu_table(norm)))>threshold], biom)
  biom <- transform_sample_counts(biom, function(x) x+pseudocounts)
  lvls <- unique(unlist(sample_data(biom)[, group]))
  diagdds <- phyloseq_to_deseq2(biom, as.formula(paste0("~",group)))
  diagdds <- DESeq(diagdds, test="Wald", fitType="mean")
  res <- results(diagdds, cooksCutoff = FALSE)
  sigtab <- res[which(res$padj < alpha), ]
  if (nrow(sigtab) > 0 ) {
    sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(biom)[rownames(sigtab), ], "matrix"))
    expr <- substitute(x != "NA", list(x = as.name(taxlevel)))
    sigtab$log2FoldChange <- round(sigtab$log2FoldChange, 3)
    sigtab$lfcSE <- round(sigtab$lfcSE, 3)
    sigtab$pvalue <- signif(sigtab$pvalue, 3)
    sigtab$padj <- signif(sigtab$padj, 3)
    if (taxlevel == "OTU") {
      sigtab$OTU <- rownames(sigtab)
    }
    data <- sigtab[,c("log2FoldChange", "lfcSE", "pvalue", "padj", taxlevel)] %>%
      arrange(padj) %>% filter_(expr)
    if (taxlevel == "OTU") {
      data$Taxon <- apply(data, 1, function(x){tax_table(biom)[x["OTU"],"Rank6"]})
    }
    data[1:min(nrow(data),maxgroups),]
  }
  else {
    data.frame()
  }
}

#' Make a boxplot of significant OTUs found with my_deseq2
#'
#' @param biom a phyloseq object.
#' @param taxlevel string giving taxonomic level (from tax_table) to plot at.
#' @param condition string giving column in sample data to compare
#' @param results tibble returned from my_deseq2
#' @param glom run tax_glom on biom to get to desired level (FALSE if already at right level)
#' @param printSig print p-values on plot
#' @param cex font size
#' @param colors optional vector of colors
#' @export
#' @examples
#' boxplot_biom()

boxplot_biom <- function(biom, taxlevel, condition, results, title=NULL,
                       glom = TRUE, printSig = TRUE, cex = 2,
                       colors = "white", show_points = TRUE, pointColors=NULL,
                       pointShapes=NULL, minPer = 10, log=FALSE) {
  library(phyloseq)
  library(ggplot2)
  nlev <- nlevels(factor(sample_data(biom)[[condition]]))
  if (is.null(colors)) {
    colors <- c("white")
  }
  group <- sample_data(biom)[[condition]]
  if (glom) {
    biom <- tax_glom(biom, taxrank = taxlevel)
  }
  norm <- transform_sample_counts(biom, function(x) 100*x / sum(x))
  if (taxlevel !="OTU") {
    taxa_names(norm) <-  make.unique(tax_table(norm)[,taxlevel])
  }
  norm <- prune_taxa(results[[taxlevel]] %>% as.character, norm)
  otus <- as.data.frame(otu_table(norm))
  otus$taxon <- row.names(otus)
  otus <- gather(otus, sample, abundance, -taxon)
  otus[[condition]] <- as.factor(sample_data(norm)[otus$sample,][[condition]])
  otus$group <- interaction(otus[[condition]],factor(otus[["taxon"]], levels=results[[taxlevel]]))
  abundant_taxa <- aggregate(abundance~taxon, otus, max) %>% .[.$abundance>minPer,]
  otus <- otus[otus$taxon %in% abundant_taxa$taxon, ]
  results <- data.frame(results)
  colors <- rep(colors, nlevels(factor(sample_data(biom)[[condition]]))*nrow(results))
  results <- results[results[,taxlevel] %in% otus$taxon,]
  if (log) {
    otus$abundance <- log2((otus$abundance*100)+1)
    ylab = "log2(Relative abundance)"
  } else {
    ylab = "Relative Abundance (%)"
  }
  p <- ggplot(otus, aes(group, abundance, fill=group))
  p <- p + geom_boxplot(outlier.colour = "white") + ylab(ylab)
  if (show_points) {
    p <- p + geom_jitter()
  }
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  xlab("") + ggtitle(title) +
  theme(panel.background = element_rect(fill = 'white', color='black'))
  if (printSig) {
    for (i in 1:nrow(results)) {
      p <- p + annotate("text", label = signif(results[i,]$padj,2) , x = nlev*i-(nlev/4.0), y= -3 , cex)
    }
  }
  p <- p + scale_fill_manual(values=colors) + guides(fill=FALSE)
  p
}

