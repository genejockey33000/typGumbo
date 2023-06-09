#' View GO Enrichment Results
#'   simple viewer of significantly enriched pathways that are output from mkEnrich()
#'
#' @param x Output object from mkEnrich()
#' @param view What to focus on for this view. One of "native" (default),
#'   "best.qvals" (most sig. qvalues)
#'   "best.CCs" (most sig. in Cellular category)
#'   "best.BPs" (most sig. in Biological Process category)
#'   "best.MFs" (most sig. in Molecular Function category)
#'
#' @return Streamlined version of results from mkEnrich object
#' @importFrom utils View
#' @export
#'
ViewEnrich <- function(x, view = "native") {
  if (!(view %in% c("native","best.qvals","best.CCs","best.BPs","best.MFs"))) {stop(
    "view should be one of ('best.qvals','best.CCs','best.BPs','best.MFs')"
  )}
  quickView <- as.data.frame(x@result[,(colnames(x@result) %in% c("ONTOLOGY","ID", "Description", "GeneRatio", "qvalue", "geneID"))])

  if (view == "native") {
    quickView <- quickView
  } else if (view == "best.qvals"){
  quickView <- quickView[order(quickView$qvalue),]
  } else if (view == "best.CCs") {
    quickView <- quickView[grepl(pattern = "CC", quickView$ONTOLOGY),]
  } else if (view == "best.BPs") {
    quickView <- quickView[grepl(pattern = "BP", quickView$ONTOLOGY),]
  } else if (view == "best.MFs") {
    quickView <- quickView[grepl(pattern = "MF", quickView$ONTOLOGY),]
  }
  View(quickView)
}
