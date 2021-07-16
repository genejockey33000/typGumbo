#' Pairwise Correlation Skew NULL Distribution (pairwiseCorSkewNULL)
#'
#'Generates a NULL distribution for "pairwiseCorSkew" to test how frequently
#' one expects to see an observed skew by chance. Note that this can take a
#' while to run if inputing large matrices with high numbers of iterations.
#' It's a good idea to test with 100 iterations first before ramping up.
#' After generation x <- pairwiseCorSkewNULL, the distribution can be visualized
#' by plotting plot(density(x)), then abline(v = observedSkew from pairwiseCorSkew)
#'
#'
#' @param x matrix 1
#' @param y matrix 2
#' @param method  Type of correlation to perform "spearman" or "pearson"
#' @param iter Number of iterations
#' @param pct percentage of upper variance to include
#'
#' @export
#'
pairwiseCorSkewNULL <- function(x, y, method = "spearman", iter = 1000, pct = 1) {
  t1 <- Sys.time()

  x <- x[,(colnames(x) %in% colnames(y))]
  x <- x[,sort(colnames(x))]
  y <- y[,(colnames(y) %in% colnames(x))]
  y <- y[,sort(colnames(y))]

  x <- x[(row.names(x) %in% row.names(y)),]
  x <- x[sort(row.names(x)),]
  y <- y[(row.names(y) %in% row.names(x)),]
  y <- y[sort(row.names(y)),]



  if (all(colnames(x) == colnames(y)) & all(row.names(x) == row.names(y))) {
    print("nice matrices!")
    if (pct < 1) {
      print("isolating specified upper variance of matrix 1 variant measures")
      bvar <- apply(x,1,stats::var)
      combo <- cbind(x,y,bvar)
      comboSort <- combo[order(-bvar),]
      comboUQ <- comboSort[c(1:(nrow(comboSort)*pct)),1:(ncol(comboSort)-1)]
      x <- t(comboUQ[,c(1:(ncol(comboUQ)/2))])
      y <- t(comboUQ[,c(((ncol(comboUQ)/2)+1):ncol(comboUQ))])
    } else {
      x<- t(x)
      y <- t(y)
    }
    cat("permuting y matrix", iter ,"times, and calculating skews.\n")
    RANDrhoSkews <- vector(mode = "numeric", length = iter)
    RANDrho <- vector(mode = "numeric", length = ncol(x))
    n <- 1
    while (n <= iter) {
      RANDy <- y[sample(nrow(y), replace = FALSE),]

      for (i in 1:(ncol(x))) {
        RANDrho[i] <- stats::cor(x[,i],RANDy[,i], method = method)
      }
      RANDrho <- sort(RANDrho, decreasing = TRUE)
      top100 <- RANDrho[1:100]
      bottom100 <- RANDrho[(length(RANDrho)-99):length(RANDrho)]
      RANDrhoSkew <- abs(sum(top100)) / abs(sum(bottom100))
      RANDrhoSkews[n] <- RANDrhoSkew
      reportProgress(n, iter, t1)
      n <- n+1
    }
    runtime <- Sys.time() - t1
    cat("It took", runtime, attr(runtime, "units"), "to run", iter, "iterations of", ncol(x), "comparisons")
    return(RANDrhoSkews)

  } else {print("Matrices are not compatable. I don't know what you did.")}
}
