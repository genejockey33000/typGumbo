#' Make Gene Ontology Enrichment Object
#' Simple script that uses clusterProfiler to generate GO enrichment object
#' Outputted object can be fed into downstream functions like gene concept network plots
#' Can take a bit of time to run
#'
#' @param x quoted path to *.csv file containing ONE column of gene names
#' @param db One of "GO" (default), or "Reactome". May add more later.
#' @param GOont Enter one of "ALL" (default), "BP", "MF", "CC"
#' @param qvalueCutoff Highest qvalue returned (default = .01)
#'
#' @return GO enrichment object
#' @importFrom clusterProfiler enrichGO
#' @export
mkEnrich <- function(x, db = "GO",GOont = "ALL", qvalueCutoff = 0.01) {
  input <- scan(file = x, what = character(), sep = ",")
  input <- unique(input)
  if (db == "Reactome") {
    input <- clusterProfiler::bitr(input, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    output <- ReactomePA::enrichPathway(input$ENTREZID, qvalueCutoff = qvalueCutoff, readable = TRUE)
  } else {
  output <- clusterProfiler::enrichGO(input, OrgDb = "org.Hs.eg.db", ont=GOont, keyType = "SYMBOL", qvalueCutoff = qvalueCutoff)
  }
  return(output)
}
