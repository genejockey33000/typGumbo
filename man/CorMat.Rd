% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CorMat.R
\name{CorMat}
\alias{CorMat}
\title{Matrix Correlations, all to all (CorMat)}
\usage{
CorMat(x, y, method = "spearman", adj.meth = "BH")
}
\arguments{
\item{x}{matrix 1}

\item{y}{matrix 2}

\item{method}{"pearson" (default) or "spearman" (rank order correlation)}

\item{adj.meth}{Method to use for multiple comparisons correction default = "BH" for Benjamini Hochberg, other options are, "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none" uses stats::p.adjust function}
}
\description{
Starting with two matrices with compatible row.names
(corresponding to *samples*)
and different column.names (corresponding to *measurements*)
i.e.
mat1 = measures of cognition
mat2 = RNAseq or Tau measurements
CorMat will sync the rows between matrices
reports surviving rows and correlates
every column in mat1 with every column in mat2
Returns a flat data frame with all comparisons,
rvals, and pvals(ranked).
}
