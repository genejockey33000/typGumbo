#' Clean limma test table
#'
#' Takes a test table outputted from limma topTable function and generates table that
#' can be easily used for ggplot functions. Can be directly input to the Vulcan function
#' for plotting volcano plots. Assumes that row.names of x are ENSG IDs. Those ENSG IDs are
#' mapped using the "map" file with column names "ensembl_gene_id" and "gene" for ENSG ID and
#' gene name respectively.
#'
#' @param x test table output from limma topTable function
#' @param lfc.cut minimum  absolute log fold change for significant hit call. Default = 0.5
#' @param adj.p.cut adjusted P value cutoff for significant hit call. Default = .05
#' @param map map that connects ENSG values (or other row.names) to common gene names
#'
#' @export
cleanTable <- function (x, lfc.cut = 0.5, adj.p.cut = 0.05, map = map) {
  #remove any rows with NA values
  chop <- apply(x, 1, function(x) {
    sum(is.na(x)) < 1
  })
  y <- x[chop, ]
  #add columns for ENSGIDs and gene names
  ensg <- row.names(y)
  ext_gene <- map$gene[match(ensg, map$ensembl_gene_id)]
  y <- cbind.data.frame(ensg, ext_gene, y)
  rownames(y) <- NULL

  #remove ENSGs without proper gene names
  chop <- y$ext_gene != ""
  y <- y[chop, ]

  #add DE column and make calls
  y$DE <- "NO"
  y$DE[y$logFC > lfc.cut & y$adj.P.Val < adj.p.cut] <- "UP"
  y$DE[y$logFC < -lfc.cut & y$adj.P.Val < adj.p.cut] <- "DOWN"

  #add rank column and sort by rank for fun
  y <- dplyr::mutate(y, rank = -log10(y$P.Value) * sign(logFC))
  y <- dplyr::arrange(y, dplyr::desc(rank))
  return(y)
}
