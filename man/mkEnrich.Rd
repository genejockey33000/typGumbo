% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mkEnrich.R
\name{mkEnrich}
\alias{mkEnrich}
\title{Make Gene Ontology Enrichment Object}
\usage{
mkEnrich(x, db = "GO", GOont = "ALL", qvalueCutoff = 0.01)
}
\arguments{
\item{x}{quoted path to *.csv file containing ONE column of gene names}

\item{db}{One of "GO" (default), or "MKEGG". May add more later but Reactome is being problematic.}

\item{GOont}{Enter one of "ALL" (default), "BP", "MF", "CC"}

\item{qvalueCutoff}{Highest qvalue returned (default = .01)}
}
\value{
GO enrichment object
}
\description{
Simple script that uses clusterProfiler to generate GO enrichment object.
Takes single column .csv file of gene symbols (i.e. "MAPT", "NEUROG2") as input (no column header).
Outputted object can be fed into downstream functions like gene concept network plots
Can take a bit of time to run.
}
