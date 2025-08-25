#' Sleuth Test table filter, removes NA rows and unmapped genes
#'
#' @param x sleuth test table from wald test
#'
#' @export
#'
noNANAtt <- function(x) {
  chop <- apply(x, 1, function(x) {sum(is.na(x)) < 1})
  y <- x[chop,]
  chop <- y$ext_gene != ""
  return(y[chop,])
}
