#' Subset vector to only EVEN indexes
#'
#' For a vector x, even(x) returns the even members. Why do you want this? There are times. Trust me.
#'
#' @param x any vector
#' @noRd
#'
even <- function(x) {
  str <- 1:(length(x))
  even <- !(as.logical(str %% 2))
  return(x[even])
}
