#' rcorr with adjusted P values included
#'
#' This is just Hmisc's rcorr function with two changes
#' 1) NA in diagonal are replaced with 0's
#' 2) returns adjusted pvalues in list object
#'
#' @param x a numeric matrix with at least 5 rows and at least 2
#'   columns (if y is absent). For print, x is an object produced
#'   by rcorr.
#' @param y a numeric vector or matrix which will be concatenated to x.
#'   If y is omitted for rcorr, x must be a matrix.
#' @param type specifies the type of correlations to compute. Spearman correlations are the Pearson linear correlations computed on the ranks of non-missing elements, using midranks for ties.
#' @param p.adj Method to be used for p value adjustment "BH" default
#'
#' @export
rcorradj <- function (x, y, type = c("pearson", "spearman"), p.adj = "BH")
{
  type <- match.arg(type)
  if (!missing(y))
    x <- cbind(x, y)
  x[is.na(x)] <- 1e+50
  storage.mode(x) <- "double"
  p <- as.integer(ncol(x))
  if (p < 1)
    stop("must have >1 column")
  n <- as.integer(nrow(x))
  if (n < 5)
    stop("must have >4 observations")
  h <- .Fortran(F_rcorr, x, n, p, itype = as.integer(1 + (type ==
                                                            "spearman")), hmatrix = double(p * p), npair = integer(p *
                                                                                                                     p), double(n), double(n), double(n), double(n), double(n),
                integer(n))
  npair <- matrix(h$npair, ncol = p)
  h <- matrix(h$hmatrix, ncol = p)
  h[h > 1e+49] <- NA
  nam <- dimnames(x)[[2]]
  dimnames(h) <- list(nam, nam)
  dimnames(npair) <- list(nam, nam)
  P <- matrix(2 * (1 - pt(abs(h) * sqrt(npair - 2)/sqrt(pmax(1 -
                                                               h * h, 0)), npair - 2)), ncol = p)
  P[abs(h) == 1] <- 0
  diag(P) <- 0
  dimnames(P) <- list(nam, nam)

  padj <- stats::p.adjust(P, method = p.adj, n = n)
  row.names(padj) <- row.names(P)
  colnames(padj) <- colnames(P)
  dim(padj) <- dim(P)
  structure(list(r = h, n = npair, P = P, P.adj = padj), class = "rcorr")
}
