#' Filter gene matrices based on k over A
#'
#' This is a slightly easier and shorter function to run the k over A algorithm.
#' Filters matrices so that output only includes measurements (genes, proteins, etc.)
#' higher than A in at least k samples.
#'
#' @param x matrix to be filtered with measurements in rows and samples in columns
#' @param k number of samples that need to be higher than A to pass filter
#' @param A minimum value to pass filter
#'
#' @return filtered matrix
#' @export
#'
chop <- function(x, k = NULL, A = NULL) {
  f1 <- genefilter::kOverA(k, A=A, na.rm = TRUE)
  flist <- genefilter::filterfun(f1)
  filt <- genefilter::genefilter(x, flist)
  out <- x[filt,]
  return(out)
}
