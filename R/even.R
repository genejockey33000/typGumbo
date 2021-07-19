#' Subset vector to only EVEN members
#'
#' For a vector x, even(x) returns the even members
#'
#' @param x any vector
#' @noRd
#'
even <- function(x) {
  str <- 1:(length(x))
  even <- !(as.logical(str %% 2))
  return(x[even])
}
