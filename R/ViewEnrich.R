#' View GO Enrichment Results
#'   simple viewer of significantly enriched pathways that are output from mkEnrich()
#'
#' @param x Output object from mkEnrich()
#'
#' @return Streamlined version of results from mkEnrich object
#' @export
#'
ViewEnrich <- function(x) {
  quickView <- as.data.frame(x@result[,c(1,3,4,8,9)])
  View(quickView)
}
