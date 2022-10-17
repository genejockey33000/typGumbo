#' Make Rank File for GSEA analysis
#'
#' Takes comparison table with ENSG IDs as row.names and columns for
#' significance and for direction
#' Ranks genes using -log10(sig) * dir
#'
#' @param x Test table with ENSGs as row.names and columns for significance
#' of test and direction of change
#' Designed to directly convert limma topTable.
#' @param fileName specify name of file to be written including .rnk suffix.
#' (i.e. file = "output.rnk")
#' @param sig name of column containing significance, not quoted
#' @param direc name of column containing directionality. b value for DEG or r value for correlation
#'
#' @return rank file
#' @export
#'
makeRank <- function(x, fileName, sig, direc) {
  rnk <- x %>%
    dplyr::mutate(rank = -log10(sig)* direc) %>%
    dplyr::mutate(ext_gene = map$ext_gene[match(row.names(x), map$ens_gene)]) %>%
    dplyr::select(ext_gene, rank) %>%
    dplyr::arrange(dplyr::desc(rank))
  utils::write.table(rnk, file = fileName, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  return(rnk)
}
