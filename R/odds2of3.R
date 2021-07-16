#' Odds of passing 2 out of 3 Comparisons
#'
#' Highly specialized function to compare whether a set of values in matrix y, are significantly correlated (uncorrected pVal)
#' to 2 out of the 3 columns of values supplied in matrix x. Runs through each column in y and correlates with each
#' of the columns in x then runs a set of tests to determine if correlation passes set rval (set to .35), that
#' rvalue signs agree, and finally that 2 out of the 3 pvals exceed the cutoff (p.cutoff). Subsequently runs a NULL distribution
#' analysis to determine how frequently the observed number of passers would be expected to occur by chance (FDR)
#'
#' @param x matrix of 3 columns to compare (could make this extensible but haven't yet)
#' @param y matrix of expression values (or proteomics or anything really)
#' @param method correlation method to use "spearman", "pearson"
#' @param iter number of iterations (can get slow if too large)
#' @param p.cutoff p value cutoff for determining significance
#' @param pct pct of upper variance to include in comparisons (default = 1 or include all)
#'
#' @return
#' @export
odds2of3 <- function(x,y, method = "pearson", iter = 100, p.cutoff = 0.05, pct = 1.0) {
  t1 <- Sys.time()
  if(all(row.names(x) == row.names(y))){  ##check that matrices are compatible

    if(pct < 1) {  # check if pct variance argument is specified
      bvar <- apply(y,2,stats::var)  # calculate variance for all columns
      y <- y[,bvar > stats::quantile(bvar, probs = (1-pct))] # subset y matrix to only measures with passing variances
    }
    pass.matrix <- NULL # establish pass.matrix (output object)
    n <- 1

    while(n <= iter) {
      repeat {
        rany <- y[sample(nrow(y)),]  # randomize y matrix (proteomics)
        if ((sum(row.names(y) == row.names(rany)))/nrow(y) < 0.25) break  # test that no more than 25% of rows match original matrix
      }
      pass <- NULL
      for (j in 1:ncol(rany)) {
        pVals <- NULL
        rVals <- NULL

        for (i in 1:ncol(x)) {
          test <- Hmisc::rcorr(x[,i], rany[,j], type = method)
          pVals <- c(pVals, test$P[1,2])
          rVals <- c(rVals, test$r[1,2])
        }
        sig.rVals <- abs(rVals) > 0.35

        if (sum(sig.rVals) < 2) {
          result <- FALSE
        } else {
          pVals <- pVals[sig.rVals]
          if (sum(pVals <= p.cutoff) < 2) {
            result <- FALSE
          } else {
            rVals <- rVals[sig.rVals]
            direction <- sign(rVals)
            if (length(rVals) == 2 && direction[1] != direction[2]) {
              result <- FALSE
            } else {
              result <- TRUE
            }
          }
        }
        pass <- c(pass, result)
      }
      pass.matrix <- rbind(pass.matrix, pass)

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

    colnames(pass.matrix) <- colnames(y)

    n.pass <- apply(pass.matrix, 2, sum)
    pct.pass <- n.pass/iter
    metrics <- rbind(n.pass, pct.pass)
    colnames(metrics) <- colnames(y)

    prot.pass <- apply(pass.matrix, 1, sum)
    prot.FDR <- (sum(prot.pass >= 19)) / iter

    t2 <- Sys.time()
    time.to.run <- t2 - t1
    return(list(pass.matrix = pass.matrix, metrics = metrics, time.to.run = time.to.run, prot.pass = prot.pass, prot.FDR = prot.FDR, p.cutoff = p.cutoff))

  } else {stop("The two matrices must be in the same order (row-wise)")}
}
