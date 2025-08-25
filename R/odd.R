#' Subset vector to only ODD indexes
#'
#' For a vector x, odd(x) returns the odd members
#'
#' @param x any vector
#' @noRd
#'
odd <- function(x) {
  str <- 1:(length(x))
  odd <- as.logical(str %% 2)
  return(x[odd])
}
