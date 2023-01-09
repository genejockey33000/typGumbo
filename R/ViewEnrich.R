#' View GO Enrichment Results
#'   simple viewer of significantly enriched pathways that are output from mkEnrich()
#'
#' @param x Output object from mkEnrich()
#' @param view What to focus on for this view. One of "best.qvals" (default),
#'   "best.CCs" (most sig. in Cellular category)
#'   "best.BPs" (most sig. in Biological Process category)
#'   "best.MFs" (most sig. in Molecular Function category)
#'
#' @return Streamlined version of results from mkEnrich object
#' @export
#'
ViewEnrich <- function(x, view = "best.qvals") {
  quickView <- as.data.frame(x@result[,c(1,3,4,8,9)])
  quickView <- quickView[order(quickView$qvalue),]
  if (view == "best.CCs") {
    quickView <- quickView[grepl(pattern = "CC", quickView$ONTOLOGY),]
  } else if (view == "best.BPs") {
    quickView <- quickView[grepl(pattern = "BP", quickView$ONTOLOGY),]
  } else if (view == "best.MFs") {
    quickView <- quickView[grepl(pattern = "MF", quickView$ONTOLOGY),]
  }
  View(quickView)
}
