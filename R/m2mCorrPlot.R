#' Dot matrix of Correlations
#'
#' Generates a matrix plot of graphical representations of correlations using either 1 or 2 matrices.
#' This is a wrapper for the corrplot function (from the CorrPlot package) that adds the step of combining
#' matrices into one (if necessary) and calculating correlations for plotting. The problem of graphing larger
#' matrices still remains. This is generally good up to 10 columns of measures after which the graphics can get
#' ugly. Generally good to minimize the length of measurement (column) names.
#'
#' @param mat1 matrix one with samples in rows and measurements in columns
#' @param mat2 optional: matrix of measures using overlapping samples as matrix 1. Samples not appearing in both
#'   matrices will be culled. mat2 will be transposed if necessary.
#' @param method Type of correlation to perform "pearson" or "spearman"
#'
#' @return
#' @export
m2mCorrPlot <- function(mat1, mat2, method = "pearson"){
  if (missing(mat2)) {
    cmat.cor <- Hmisc::rcorr(mat1, type = method)
    cmat.cor$P[is.na(cmat.cor$P)] <- 0
  } else {
    cmat <- stickem(mat1, mat2)
    cmat.cor <- Hmisc::rcorr(cmat, type = method)
    cmat.cor$P[is.na(cmat.cor$P)] <- 0
  }
  return(corrplot::corrplot(cmat.cor$r, method = "circle", order = 'original', type = "upper", diag = FALSE,p.mat = cmat.cor$P, sig.level = 0.05, insig = 'label_sig'))
}
