#' Version Out for removing ENSG and ENST version extensions
#'
#' In some gene or transcript expression data there are trailing numbers
#'     after a decimal indicating the version of that sequence originally
#'     matched against. E.g. ENSG000000876.5  - these can cause problems
#'     with mapping to gene meta data and are often best removed for downstream
#'     analysis
#'
#' @param x a vector of transcripts(ENSTs) or genes(ENSGs) with versions
#'
#' @export
#'
verOut <- function(x) {
  return(sub("\\.[0-9]*", "", x))
}
