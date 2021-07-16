#' Pairwise correlation matrix
#'
#' Temporary description
#'
#' @param x matrix 1
#' @param y matrix 2
#' @param method "pearson" or "spearman"
#' @param pct percentage of upper variance to include
#'
#' @export
pairwiseCorMat <- function(x,y, method = "pearson", pct = "1.0")  {
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
    if (pct < 1.0) {
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
    cat("Calculating correlations for", nrow(x), "measurements across", ncol(x), "samples")
    CWC <- NULL
    for (i in 1:ncol(x))  {
      xname <- character()
      rval <- numeric()
      pval <- numeric()
      v1 <- x[,i]
      v2 <- y[,i]
      output <- Hmisc::rcorr(v1, v2, type = method)
      xname <- colnames(x)[i]
      rval <- c(rval, output$r[1,2])
      pval <- c(pval, output$P[1,2])

      Cors <- cbind.data.frame(xname, rval, pval)
      CWC <- rbind.data.frame(CWC,Cors)
    }
    CWC <- CWC[order(CWC$pval),]
    skew <- sum(CWC$rval)
    CWC <- list("CWC" = CWC, "skew" = skew)
    return(CWC)
  } else {print("Matrices are incompatable. I don't know what you did.")}
}
