#' Run DESeq2 on a phyloseq object and return tibble of significantly abundant OTUs
#'
#' @param biom a phyloseq object.
#' @param taxlevel string giving taxonomic level (from tax_table) to compare at.
#' @param condition string giving column in sample data to compare
#' @param alpha p-value cutoff
#' @param threshold minimum percentage of OTU to report
#' @param glom run tax_glom on biom to get to desired level (FALSE if already at right level)
#' @export
#' @examples
#' my_deseq2()

my_deseq2 <- function(biom, taxlevel, condition, alpha = 0.01, threshold=1, glom=TRUE) {
  if (glom) {
    biom <- tax_glom(biom, taxrank = taxlevel)
  }
  norm <-  transform_sample_counts(biom, function(x) 100*x / sum(x))
  biom <- prune_taxa(taxa_names(norm)[rowMeans(as.data.frame(otu_table(norm)))>threshold], biom)
  lvls <- levels(unlist(sample_data(biom)[, condition]))
  diagdds <- phyloseq_to_deseq2(biom, as.formula(paste0("~",condition)))
  diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")
  res <- results(diagdds, cooksCutoff = FALSE)
  sigtab <- res[which(res$padj < alpha), ]
  sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(biom)[rownames(sigtab), ], "matrix"))
  sigtab$direction <- ifelse(sigtab$log2FoldChange > 0,  lvls[2], lvls[1])
  expr <- substitute(x != "NA", list(x = as.name(taxlevel)))
  sigtab$log2FoldChange <- round(sigtab$log2FoldChange, 3)
  sigtab$lfcSE <- round(sigtab$lfcSE, 3)
  sigtab$pvalue <- signif(sigtab$pvalue, 3)
  sigtab$padj <- signif(sigtab$padj, 3)
  select(as.tibble(sigtab), c("log2FoldChange", "lfcSE", "pvalue", "padj", taxlevel,
                              "direction")) %>% arrange(padj) %>% filter_(expr)
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
#' my_boxplot()

my_boxplot <- function(biom, taxlevel, condition, results, title=NULL,
                       glom = TRUE, printSig = TRUE, cex = 2, colors = rep("white",2*nrow(results))) {
  group <- sample_data(biom)[[condition]]
  if (glom) {
    biom <- tax_glom(biom, taxrank = taxlevel)
  }
  norm <- transform_sample_counts(biom, function(x) 100*x / sum(x))
  taxa_names(norm) <- make.unique(as.character(as.data.frame(tax_table(norm))[[taxlevel]]))
  norm <- prune_taxa(results[[taxlevel]] %>% as.character, norm)
  otus <- as.data.frame(otu_table(norm))
  otus$taxon <- row.names(otus)
  otus <- gather(otus, sample, abundance, -taxon)
  otus[[condition]] <- as.factor(sample_data(norm)[otus$sample,][[condition]])
  otus$group <- interaction(otus[[condition]],factor(otus[["taxon"]], levels=results[[taxlevel]]))
  p = ggplot(otus, aes(group, abundance, fill=group)) +
    geom_boxplot(outlier.colour = "white") + ylab("Relative Abundance (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ggtitle(title) +
    theme(panel.background = element_rect(fill = 'white', color='black'))
  if (printSig) {
    for (i in 1:nrow(results)) {
      p <- p + annotate("text", label = signif(results[i,]$padj,2) , x = 2*i-0.5, y= -3 , cex)
    }
  }
  p <- p + scale_fill_manual(values=colors) + guides(fill=FALSE)
  p
}
