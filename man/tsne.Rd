% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsne.R
\name{tsne}
\alias{tsne}
\title{TSNE Plot Maker: Copied from M3C}
\usage{
tsne(
  mydata,
  labels = FALSE,
  perplex = 15,
  printres = FALSE,
  seed = FALSE,
  axistextsize = 18,
  legendtextsize = 18,
  dotsize = 5,
  textlabelsize = 4,
  legendtitle = "Group",
  controlscale = FALSE,
  scale = 1,
  low = "grey",
  high = "red",
  colvec = c("skyblue", "gold", "violet", "darkorchid", "slateblue", "forestgreen",
    "violetred", "orange", "midnightblue", "grey31", "black"),
  printheight = 20,
  printwidth = 22,
  text = FALSE
)
}
\arguments{
\item{mydata}{Data frame or matrix: if dataframe/matrix should have samples as columns and rows as features}

\item{labels}{Character vector: if we want to just label by sex for example}

\item{perplex}{Numerical value: perplexity value that Rtsne uses internally}

\item{printres}{Logical flag: whether to print the t-SNE into current directory}

\item{seed}{Numerical value: optionally set the seed}

\item{axistextsize}{Numerical value: axis text size}

\item{legendtextsize}{Numerical value: legend text size}

\item{dotsize}{Numerical value: dot size}

\item{textlabelsize}{Numerical value: text inside plot label size}

\item{legendtitle}{Character vector: text legend title}

\item{controlscale}{Logical flag: whether to control the colour scale}

\item{scale}{Numerical value: 1=spectral palette, 2=manual low and high palette, 3=categorical labels}

\item{low}{Character vector: continuous scale low color}

\item{high}{Character vector: continuous scale high color}

\item{colvec}{Character vector: a series of colors in vector for categorical labels, e.g. c("sky blue", "gold")}

\item{printheight}{Numerical value: png height}

\item{printwidth}{Numerical value: png width}

\item{text}{Character vector: if we wanted to label the samples with text IDs to look for outliers}
}
\description{
This is just a copy of the very nice tsne wrapper from the M3C package which uses the Rtsne
This is a flexible t-SNE function that can be run on a standard data frame.
It is a wrapper for Rtsne/ggplot2 code and can be customized with different
colors and font sizes and more.
}
