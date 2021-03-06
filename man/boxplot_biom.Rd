% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deseq2_biom.R
\name{boxplot_biom}
\alias{boxplot_biom}
\title{Make a boxplot of significant OTUs found with my_deseq2}
\usage{
boxplot_biom(biom, taxlevel, condition, results, title = NULL,
  glom = TRUE, printSig = TRUE, cex = 2, colors = "white",
  show_points = TRUE, pointColors = NULL, pointShapes = NULL,
  minPer = 10, log = FALSE)
}
\arguments{
\item{biom}{a phyloseq object.}

\item{taxlevel}{string giving taxonomic level (from tax_table) to plot at.}

\item{condition}{string giving column in sample data to compare}

\item{results}{tibble returned from my_deseq2}

\item{glom}{run tax_glom on biom to get to desired level (FALSE if already at right level)}

\item{printSig}{print p-values on plot}

\item{cex}{font size}

\item{colors}{optional vector of colors}
}
\description{
Make a boxplot of significant OTUs found with my_deseq2
}
\examples{
boxplot_biom()
}
