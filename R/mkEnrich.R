#' Make Gene Ontology Enrichment Object
#' Simple script that uses clusterProfiler to generate GO enrichment object
#' Outputted object can be fed into downstream functions like gene concept network plots
#' Can take a bit of time to run
#'
#' @param x quoted path to *.csv file containing ONE column of gene names
#' @param ont Enter one of "ALL" (default), "BP", "MF", "CC"
#' @param qvalueCutoff Highest qvalue returned (default = .01)
#'
#' @return GO enrichment object
#' @importFrom clusterProfiler enrichGO
#' @export
mkEnrich <- function(x, ont = "ALL", qvalueCutoff = 0.01) {
  input <- scan(file = x, what = character(), sep = ",")
  inputGO <- clusterProfiler::enrichGO(input, OrgDb = "org.Hs.eg.db", ont=ont, keyType = "SYMBOL", qvalueCutoff = qvalueCutoff)
  return(inputGO)
}
