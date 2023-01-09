#' View GO Enrichment Results
#'   simple viewer of significantly enriched pathways that are output from mkEnrich()
#'
#' @param x Output object from mkEnrich()
#'
#' @export
#'
ViewEnrich <- function(x) {
  out <- as.data.frame(x@result)
  View(out[,c(1,3,4,8,9)])
}
