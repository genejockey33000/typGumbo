#' Subset vector to only ODD members
#'
#' For a vector x, odd(x) returns the odd members
#'
#' @param x any vector
#'
#' @export
odd <- function(x) {
  str <- 1:(length(x))
  odd <- as.logical(str %% 2)
  return(x[odd])
}
