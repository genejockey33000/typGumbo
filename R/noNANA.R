#' Removes rows with any NA values
#'
#' This is just a quicker function to remove any row containing any
#' NA value. Note this may be overkill in many cases. I use it for
#' processing proteomic matrices where NA is an indicator of unreliable
#' measurements.
#'
#' @param x matrix with NA values
#'
#' @return
#' @export
#'
NoNANA <- function(x) {
  chop <- apply(x, 1, function(x) {sum(is.na(x)) < 1})
  return(x[chop,])
}
