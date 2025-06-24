#' Clean Sleuth Test Table
#'
#' Cleans up test table output from sleuth_results.
#'
#' @param x Test table from sleuth_results
#' @param bcut Cutoff for fold change (z score) significance
#' @param qcut Cutoff for qvalue (FDR) significance
#'
#' @return cleaned test table
#' @export
cleanTT <- function(x, bcut = .5, qcut = .05, TT.type = "limma") {
  colnames(x)[grepl("ext_gene")] <- "gene"
  colnames(x)[grepl("^qval$")] <- "adj.P.Val"
  colnames(x)[grepl("^b$")] <- "logFC"
  chop <- apply(x, 1, function(x) {sum(is.na(x)) < 1})
  y <- x[chop,]
  chop <- y$ext_gene != ""
  d <- y[chop,]
  d$DE <- "NO"
  d$DE[d$b > bcut & d$qval < qcut] <- "UP"
  d$DE[d$b < -bcut & d$qval < qcut] <- "DOWN"
  d <- dplyr::mutate(d, rank = -log10(x$qval)*x$b)
  d <- dplyr::arrange(d, dplyr::desc(rank))
  return(d)
}
