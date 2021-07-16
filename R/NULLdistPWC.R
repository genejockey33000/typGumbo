#' The Null Distribution generator for pairwiseCorMat
#'
#' Generates a Null distribution and tests significance of skew seen in pairwiseCorMat
#'
#' @param x matrix 1
#' @param y matrix 2
#' @param method c("spearman", "pearson")
#' @param iter number of iterations (over 1000 and this gets very slow)
#' @param pct percentage of upper variance genes to test
#'
#' @export
#'
NULLdistPWC <- function(x, y, method = "spearman", iter = 1000, pct = 1) {
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
    cat("permuting y matrix", iter ,"times, and calculating skews")
    RANDrho <- NULL
    RANDrhoSUM <- NULL
    RANDrhoSUMs <- NULL

    cn <- ncol(x)
    n <- 1
    while (n<=iter) {
      RANDy <- y[sample(nrow(y), replace = FALSE),]
      ccn <- 1

      while (ccn<=cn) {
        RANDrho <- rbind(RANDrho, stats::cor(x[,ccn],RANDy[,ccn], method = method))
        ccn <- ccn+1
      }
      RANDrho <- sort(RANDrho, decreasing = TRUE)
      RANDrhoSUM <- sum(RANDrho)
      RANDrhoSUMs <- rbind(RANDrhoSUMs, RANDrhoSUM)
      RANDrho <- NULL

      if(n == .25 * iter) {
        t.25 <- (Sys.time() - t1)
        cat(t.25, attr(t.25, "units"),"have elapsed: 25% iterations complete","\n")}
      if(n == .5 * iter) {
        t.5 <- (Sys.time() - t1)
        cat(t.5, attr(t.5, "units"), "have elapsed: 50% iterations complete","\n")}
      if(n == .75 * iter) {
        t.75 <- (Sys.time() - t1)
        cat(t.75, attr(t.75, "units"), "have elapsed: 75% iterations complete","\n")}

      n <- n+1
    }
    return(RANDrhoSUMs)

  } else {print("Matrices are not compatable. I don't know what you did.")}
}
