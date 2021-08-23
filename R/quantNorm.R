#' Quantile Normalization of Expression (RNAseq) data
#'
#' Just a simple wrapper for preprocessCore() which erases colnames and row.names
#' during execution
#'
#' @param x Matrix to be quantile normalized
#'
#' @return
#' @export
#'
quantNorm <- function(x) {
  rns <- row.names(x)
  samps <- colnames(x)
  outq <- preprocessCore::normalize.quantiles(x)
  row.names(outq) <- rns
  colnames(outq) <- samps
  return(outq)
}
