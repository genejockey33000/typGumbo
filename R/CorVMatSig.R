#' Correlate vector to matrix and test significance
#'
#' Correlate Vector to rows in a  matrix and determine significance against a NULL distribution. Input vector (x)
#'  is a named vector that corresponds to row.names in matrix (y). Correlations between
#'
#' @param x Named Vector to be correlated with matrix y
#' @param y Matrix with compatible names as vector x
#' @param method type of correlation to use "pearson" or "spearman"
#' @param iter number of iterations to run
#' @param p.cutoff p value cutoff for determining significance
#' @param pct percentage upper variance to use
#'
#' @export
#'
CorVMatSig <- function(x,y, method = "pearson", iter = 100, p.cutoff = 0.05, pct = 1.0) {
  t1 <- Sys.time()
  x <- x[names(x) %in% row.names(y)]
  y <- y[row.names(y) %in% names(x),]
  x <- x[match(row.names(y), names(x))]
  if(length(x) < 6){stop("The need to be 6 or more common measurements between
                         your vector and your matix to produce meaningful correlations")}
  if(all(names(x) == row.names(y))){
    if(pct < 1) {
      bvar <- apply(y,2,stats::var)
      y <- y[,bvar > stats::quantile(bvar, probs = (1 - pct))]
    }
    rval <- NULL
    pval <- NULL
    for (j in 1:ncol(y))  {
      output <- Hmisc::rcorr(x, y[,j], type = method)
      pval <- c(pval, output$P[1,2])
      rval <- c(rval, output$r[1,2])
    }
    chop <- pval <= p.cutoff
    numSigTest <- sum(chop)

    sig.pvals <- pval[chop]
    sig.rvals <- rval[chop]
    passers <- colnames(y)[chop]
    passers <- cbind.data.frame(passers,sig.pvals, sig.rvals)
    passers <- passers[order(passers[,2]),]
    cat(numSigTest,"comparisons passed criteria in original data","\n")
    prot.sig.NULL <- NULL
    n <- 1
    print("Generating NULL Distribution")
    while(n <= iter) {
      repeat {
        ranx <- sample(x, length(x), replace = FALSE)
        if ((sum(names(x) == names(ranx)))/length(x) < 0.02) break
      }

      RANDp <- NULL
      for (j in 1:ncol(y)) {
        RANDout <- Hmisc::rcorr(ranx, y[,j], type = method)
        RANDp <- c(RANDp, RANDout$P[1,2])
      }

      protNULL <- RANDp <= p.cutoff
      prot.sig.NULL <- rbind(prot.sig.NULL, protNULL)

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
    colnames(prot.sig.NULL) <- colnames(y)
    numSigNULL <- apply(prot.sig.NULL, 1, sum)

    pct.Chance <- (sum(numSigNULL >= numSigTest)) / (iter)
    FDR <- paste("Rate by chance =", pct.Chance)
    print(FDR)
    runtime <- Sys.time() - t1
    cat("It took", runtime, attr(runtime, "units"), "to run", iter, "iterations of", ncol(y), "comparisons")
    return(list(pval = pval, numSigTest = numSigTest, numSigNULL = numSigNULL, prot.sig.NULL = prot.sig.NULL ,pct.Chance = pct.Chance, passers = passers))

  } else {stop("The named vector must match matrix row names")}

}
