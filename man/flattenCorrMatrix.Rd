% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flattenCorrMatrix.R
\name{flattenCorrMatrix}
\alias{flattenCorrMatrix}
\title{Flat matrix from Hmisc::rcorr list output
Takes the r value and the pvalue objects for the rcorr list object
  and returns a flat matrix for all combinations. Then sorts for
  best pvals}
\usage{
flattenCorrMatrix(rcorrOUT)
}
\arguments{
\item{rcorrOUT}{Correlation matrix output from rcorr or rcorr adjust}
}
\value{
flat matrix
}
\description{
Flat matrix from Hmisc::rcorr list output
Takes the r value and the pvalue objects for the rcorr list object
  and returns a flat matrix for all combinations. Then sorts for
  best pvals
}
