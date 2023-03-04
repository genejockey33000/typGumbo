#' Make Gene Ontology Enrichment Object
#' Simple script that uses clusterProfiler to generate GO enrichment object.
#' Takes single column .csv file of gene symbols (i.e. "MAPT", "NEUROG2") as input (no column header).
#' Outputted object can be fed into downstream functions like gene concept network plots
#' Can take a bit of time to run.
#'
#' @param csv quoted path to *.csv file containing ONE column of gene names
#' @param db One of "GO" (default), or "Reactome". May add more later.
#' @param GOont Enter one of "ALL" (default), "BP", "MF", "CC"
#' @param qvalueCutoff Highest qvalue returned (default = .01)
#'
#' @return GO enrichment object
#' @importFrom clusterProfiler bitr
#' @importFrom clusterProfiler enrichGO
#' @importFrom ReactomePA enrichPathway
#' @importFrom Matrix index
#' @export
mkEnrich <- function(csv, db = "GO", GOont = "ALL", qvalueCutoff = 0.01) {

  input <- scan(file = csv, what = character(), sep = ",")
  input <- unique(input)

  if (db == "Reactome") {
    input <- clusterProfiler::bitr(input, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    output <- ReactomePA::enrichPathway(input$ENTREZID, qvalueCutoff = qvalueCutoff, readable = TRUE)
  }
  if (db == "GO") {
    output <- clusterProfiler::enrichGO(input, OrgDb = "org.Hs.eg.db", ont=GOont, keyType = "SYMBOL", qvalueCutoff = qvalueCutoff)
  }
  return(output)
}
