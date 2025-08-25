#' Make a linear regression dot plot of two vectors in a matrix
#'
#' @param mat Single matrix input
#' @param x quoted name of column data for x axis (i.e. "Abeta42_40)
#' @param y quoted name of column data for y axis (i.e. "Tau-p")
#' @param datatype quoted. What kind of data is this? (i.e. "RNAseq", "Proteomics", "Assay")
#' @param title quoted. Title of plot to generate
#' @param method Statistical method for correlation reporting ("pearson" or "spearman")
#'
#' @return regression plot with regression line
#' @export
corrdot <- function(mat, x, y, datatype, title, method = "pearson") {
  x1 <- as.numeric(mat[,x])
  y1 <- as.numeric(mat[,y])
  xycor <- Hmisc::rcorr(x1,y1, type = method)
  p <- plot(x1,y1, xlab = paste0(x, " ", datatype , "; P = ",signif(xycor$P[1,2], digits = 3) ), ylab = paste0(y, " ",datatype, "; r = ", signif(xycor$r[1,2], digits = 3)), main = title)
  return(p +
           graphics::abline(stats::lm(y1~x1), col = "black"))
}
