#' Expression matrix cleaning tool
#'
#' This is a short pipeline that filters expression matrices based on k over A
#' arguments, log2() transforms and quantile normalizes. Unnecessary for ratiometric
#' data or previously processed data but useful for raw or TPM matrices.
#' Function will remove all measurements with 0 values
#' in all samples (zero variance). It will test matrix for asymmetric negative values (often resulting from
#' regressing batch effects) and raise all values by abs(smallest negative value) to remove any negative
#' expression measurement.
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
  if (min(x) < 0 & max(x) / (abs(min(x))) > 100) {
    #for matrices with unintended negative values after regression.
    #Lifts all values by abs() of lowest (negative) value
  x <- x + base::abs(min(x))
  }
  if(k > 0) {
    x <- filt(x, k = k, A = A)
    if(nrow(x) == 0) stop("Those filtering parameters removed all measurements.
                            It's doubtful you want this.")
  }
  x <- log2(x + 0.1)
  x <- quantNorm(x)
  return(x)
}
