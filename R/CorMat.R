#' @title Matrix Correlations, all to all (CorMat)
#'
#' @description Starting with two matrices with compatible row.names
#' (corresponding to *samples*)
#' and different column.names (corresponding to *measurements*)
#' i.e.
#' mat1 = measures of cognition
#' mat2 = RNAseq or Tau measurements
#' CorMat will sync the rows between matrices
#' reports surviving rows and correlates
#' every column in mat1 with every column in mat2
#' Returns a flat data frame with all comparisons,
#' rvals, and pvals(ranked).
#'
#' @param x matrix 1
#' @param y matrix 2
#' @param method "pearson" (default) or "spearman" (rank order correlation)
#' @param adj.meth Method to use for multiple comparisons correction default = "BH" for Benjamini Hochberg, other options are, "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none" uses stats::p.adjust function
#'
#' @importFrom Hmisc rcorr
#' @export
CorMat <- function(x, y, method = "spearman", adj.meth = "BH")  {
  if (is.null(dim(x))) {
    if (length(x) != nrow(y)) stop("The Vector x be the same length
      as each column of data in y. Recheck your input vector and your comparison matrix and
      try again.")
    cat("Comparing vector with",length(x)," samples to", ncol(y), "additional measurements")
    CWC <- NULL
    yname <- character()
    rval <- numeric()
    pval <- numeric()
    p.adj <- numeric()
      for (j in 1:ncol(y))  {
        v2 <- y[,j]
        ytemp <- colnames(y)[j]
        output <- Hmisc::rcorr(x, v2, type = method)
        yname <- c(yname, ytemp)
        rval <- c(rval, output$r[1,2])
        pval <- c(pval, output$P[1,2])
        p.adj <- stats::p.adjust(pval, method = adj.meth, n = ncol(y))
      }
      Cors <- cbind.data.frame(yname, rval, pval, p.adj, adj.meth)
      CWC <- rbind.data.frame(CWC,Cors)
    CWC <- CWC[order(CWC$p.adj),]
    return(CWC)
  } else {
  chop <- row.names(x) %in% row.names(y)
  x <- x[chop,]
  chop <- row.names(y) %in% row.names(x)
  y <- y[chop,]
  x <- x[sort(row.names(x)),]
  y <- y[sort(row.names(y)),]
  if (all(row.names(x) == row.names(y))) {
    print("nice matrices!")
    cat(nrow(x), "samples will be compared across", ncol(y), "measurements", ncol(x), " times!")
    CWC <- NULL
    for (i in 1:ncol(x))  {
      yname <- character()
      rval <- numeric()
      pval <- numeric()
      p.adj <- numeric()

      xtemp <- colnames(x)[i]
      v1 <- x[,i]

      for (j in 1:ncol(y))  {
        v2 <- y[,j]
        ytemp <- colnames(y)[j]
        output <- Hmisc::rcorr(v1, v2, type = method)
        yname <- c(yname, ytemp)
        rval <- c(rval, output$r[1,2])
        pval <- c(pval, output$P[1,2])
        p.adj <- stats::p.adjust(pval, method = adj.meth, n = ncol(x)*ncol(y))
      }
      xname <- character(ncol(y))
      xname[] <- xtemp
      Cors <- cbind.data.frame(xname, yname, rval, pval, p.adj, adj.meth)
      CWC <- rbind.data.frame(CWC,Cors)
      #CWC[[colnames(x)[i]]] <- Cors
    }
    CWC <- CWC[order(CWC$p.adj),]
    return(CWC)
  } else {print("Matrices are incompatable. I don't know what you did.")}
  }
}
