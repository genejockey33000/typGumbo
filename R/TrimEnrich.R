#' Trim mkEnrich() results
#'
#' Many of the enriched pathways from GO enrichment are highly overlapping.
#' i.e. BC: Ribosome and BP Ribosome Function, etc.
#' To clean up the downstream Gene concept networks it can be helpful to trim
#' redundant enriched pathways.
#'
#' @param x Output object from mkEnrich() function
#' @param terms Terms in pathways to keep. Can use multiple terms separated by pipes.
#'    i.e. "synapse|amyloid|insulin" will subset to include only terms with one of those
#'    names in the description. Not sensitive to case.
#' @return Modified enrichResult only including indicated terms
#' @export
#'
TrimEnrich <- function(x, terms=NULL) {
  if (is.null(terms)){stop("Please specify the terms you want to trim to m'kay?")}
  result <- x@result[x@result$Description[grepl(terms, x@result$Description, ignore.case = TRUE)],]
  x@result <- result
  return(x)
}
