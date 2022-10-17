#' Flat matrix from Hmisc::rcorr list output
#' Takes the r value and the pvalue objects for the rcorr list object
#'   and returns a flat matrix for all combinations. Then sorts for
#'   best pvals
#'
#' @param rcorrOUT Correlation matrix output from rcorr or rcorr adjust
#'
#' @return flat matrix
#' @export
flattenCorrMatrix <- function(rcorrOUT) {
  cormat <- rcorrOUT$r
  pmat <- rcorrOUT$P
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
