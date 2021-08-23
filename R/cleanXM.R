#' Expression matrix cleaning tool
#'
#' This is a short pipeline that filters expression matrices based on k over A
#' arguments, log2() transforms and quantile normalizes
#'
#' @param x Input expression matrix, measurements in rows
#' @param k Number of samples whose expression need to be higher than A
#' @param A Expression (measurement) level to target
#'
#' @return
#' @export
#'
cleanXM <- function(x, k = 0, A = 0) {
  x <- x[(apply(x,1,var) > 0),] #eliminates measurements with all zeros
  x <- x + base::abs(min(x)) #for matrices with negative values. Lifts all values by lowest value
  if(k > 0) {
    x <- filt(x, k = k, A = A)
    if(nrow(x) == 0) stop("Those filtering parameters removed all measurements.
                            It's doubtful you want this.")
  }
  x <- log2(x + 0.1)
  x <- quantNorm(x)
  return(x)
}
