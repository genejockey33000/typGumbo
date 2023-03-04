#' @title Make Gene Ontology Enrichment Object
#'
#' @description Simple script that uses clusterProfiler to generate GO enrichment object.
#' Takes single column .csv file of gene symbols (i.e. "MAPT", "NEUROG2") as input (no column header).
#' Outputted object can be fed into downstream functions like gene concept network plots
#' Can take a bit of time to run
#'
#' @param x quoted path to *.csv file containing ONE column of gene names
#' @param db One of "GO" (default), "Reactome", "KEGG", or "MKEGG". May add more later.
#' @param GOont Enter one of "ALL" (default), "BP", "MF", "CC"
#' @param qvalueCutoff Highest qvalue returned (default = .01)
#'
#' @return GO enrichment object
#' @importFrom clusterProfiler enrichGO
#' @export
mkEnrich <- function(x, db = "GO",GOont = "ALL", qvalueCutoff = 0.01) {
  if (db %in% c("Reactome", "KEGG", "MKEGG", "GO")) {
    cat("Preparing ", db, " analysis!")
  } else {stop("Srrrry! db argument must be either 'GO' (default), 'Reactome', 'KEGG', or 'MKEGG' (Module KEGG).")}
  input <- scan(file = x, what = character(), sep = ",")
  input <- unique(input)
  if (db == "Reactome") {
    input <- clusterProfiler::bitr(input, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    output <- ReactomePA::enrichPathway(input$ENTREZID, qvalueCutoff = qvalueCutoff, readable = TRUE)
  } else if (db == "KEGG") {
    input <- clusterProfiler::bitr(input, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
    input <- clusterProfiler::bitr_kegg(input$ENTREZID, fromType = "ncbi-geneid", toType = "kegg", organism = "hsa")
    output <- clusterProfiler::enrichKEGG(input$kegg, organism = "hsa", keyType = "kegg", pAdjustMethod = "BH", qvalueCutoff = qvalueCutoff)
  } else if (db == "MKEGG") {
    input <- clusterProfiler::bitr(input, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
    input <- clusterProfiler::bitr_kegg(input$ENTREZID, fromType = "ncbi-geneid", toType = "kegg", organism = "hsa")
    output <- clusterProfiler::enrichMKEGG(input$kegg, organism = "hsa", keyType = "kegg", pAdjustMethod = "BH", qvalueCutoff = qvalueCutoff)
  } else if (db == "GO") {
  output <- clusterProfiler::enrichGO(input, OrgDb = "org.Hs.eg.db", ont=GOont, keyType = "SYMBOL", qvalueCutoff = qvalueCutoff)
  }
  return(output)
}
