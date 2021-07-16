#' Comprehensive test of Pair-wise Correlation Parameters
#'
#' Run pairwiseCorMat.R function varying all parameters to
#' test the effect on skew
#'
#' @param x matrix 1
#' @param y matrix 2 (see pairwise CorMat for guidence)
#'
#' @return Data frame of skews for different adjustments
#' @export
#'
CompPairwiseCorMat <- function(x, y) {
  output <- data.frame()
  for (m in c("pearson", "spearman")) {
    for (p in c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)) {
      run <- paste(m,"-",p, sep = "")
      result <- pairwiseCorMat(x, y, method = m, pct = p)
      skew <- c(run, result$skew)
      output <- rbind.data.frame(output, skew)
    }
  }
  colnames(output) <- c("test", "skew")
  return(output)
}
