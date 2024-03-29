% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrdot.R
\name{corrdot}
\alias{corrdot}
\title{Make a linear regression dot plot of two vectors in a matrix}
\usage{
corrdot(mat, x, y, datatype, title, method = "pearson")
}
\arguments{
\item{mat}{Single matrix input}

\item{x}{quoted name of column data for x axis (i.e. "Abeta42_40)}

\item{y}{quoted name of column data for y axis (i.e. "Tau-p")}

\item{datatype}{quoted. What kind of data is this? (i.e. "RNAseq", "Proteomics", "Assay")}

\item{title}{quoted. Title of plot to generate}

\item{method}{Statistical method for correlation reporting ("pearson" or "spearman")}
}
\value{
regression plot with regression line
}
\description{
Make a linear regression dot plot of two vectors in a matrix
}
