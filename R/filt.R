#' Filter Data
#'
#' Just a wrapper for k over A filtering with k being the number of samples
#' with values over A required to not be filtered out.
#'
#' @param x Input expression matrix
#' @param k Number of samples required to have values over that specified in A
#' @param A Target value above which is required in k samples to be retained
#'
#' @return
#' @export
#'
filt <- function(x, k = NULL, A = NULL) {
  f1 <- genefilter::kOverA(k, A=A, na.rm = TRUE)
  flist <- genefilter::filterfun(f1)
  filt <- genefilter::genefilter(x, flist)
  out <- x[filt,]
  return(out)
}
