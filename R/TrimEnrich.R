#' Trim mkEnrich() results
#'
#' Many of the enriched pathways from GO enrichment are highly overlapping.
#' i.e. BC: Ribosome and BP Ribosome Function, etc.
#' To clean up the downstream Gene concept networks it can be helpful to trim
#' redundant enriched pathways.
#'
#' @param x Output object from mkEnrich() function
#' @param modules Row numbers (indices) of pathways to keep
#'
#' @return Modified enrichResult only including indicated modules
#' @export
#'
TrimEnrich <- function(x, modules=NULL) {
  if (is.null(modules)){stop("Please specify the indices (row numbers) of the modules you want to trim to m'kay?")}
  keepers <- paste0("GO:",  modules)
  result <- x@result[row.names(x@result) %in% keepers,]
  x@result <- result
  return(x)
}
