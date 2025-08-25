#' Removes rows with any NA values
#'
#' This is just a quicker function to remove any row containing any
#' NA value. Note this may be overkill in many cases. I use it for
#' processing proteomic matrices where NA is an indicator of unreliable
#' measurements.
#'
#' @param x matrix with NA values
#' @param margin indicator of what should be done when NA is encountered.
#'  "rows" (default) = delete entire row with NA
#'  "columns" = delete entire column with NA
#'  "cellempty" = replace NA value with empty space
#'  "cellzero" = replace NA value with numeric zero value
#'
#' @export
#'
noNANA <- function(x, margin = "rows") {
  if (margin == "columns") {
    chop <- apply(x, 2, function(x) {sum(is.na(x)) < 1})
    return(x[,chop])
  } else if (margin == "rows"){
  chop <- apply(x, 1, function(x) {sum(is.na(x)) < 1})
  return(x[chop,])
  } else if (margin == "cellempty") {
    x[is.na(x)] <- ""
    return(x)
  } else if (margin == "cellzero") {
    x[is.na(x)] <- 0
    return(x)
  } else stop('margin needs to be set to either "rows" (default), "columns", "cellempty", or "cellzero".\n I\'m not a miracle worker')
}
