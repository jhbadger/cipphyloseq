% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bar2.R
\name{plot_bar2}
\alias{plot_bar2}
\title{More flexible version of phyloseq's plot_bar}
\usage{
plot_bar2(biom, fill = NULL, title = NULL, colors = NULL,
  minPercent = NULL, pickTaxa = NULL, excludeTaxa = NULL,
  border_color = "black", log = FALSE, ppm = FALSE, rbiom = FALSE,
  showPercents = NULL, nomisc = FALSE)
}
\arguments{
\item{biom}{a phyloseq object.}

\item{fill}{string giving taxonomic level (from tax_table) to plot at.}

\item{title}{title of plot}

\item{colors}{vector of color names to use}

\item{minPercent}{optional limit to exclude low abundant taxa}

\item{pickTaxa}{optional list of taxa to limit graph to.}

\item{excludeTaxa}{optional list of taxa to exclude}

\item{border_color}{border color of bars}

\item{log}{log transform data}

\item{ppm}{translate data into parts per million rather than percent}

\item{rbiom}{return a phyloseq data object}

\item{showPercents}{limit graph to c(min,max) percent}
}
\description{
More flexible version of phyloseq's plot_bar
}
\examples{
plot_bar2()
}
