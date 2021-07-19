#' pairwiseCorSkew
#'
#' Takes two measurement matrices (generally expression matrices) containing measurements from sets of
#' samples to be check for congruence. For example, to compare transcriptomes from stem cell derived
#' cell cultures to transcriptomes from the same human subjects. In this case each column is a human subject
#' and each row a gene expression value. Column names must use the same IDs between matix 1 and 2 but need
#' not be in the same order. Some missing samples and measurements are allowed (the script will remove them).
#' pairwiseCorSkew returns a list with 2 components: 1) CWC is the column-wise correlation output exported as
#' a data frame containing all measurements rvals, pvals, fdr (Benjamini-Hochberg), and family-wise error
#' rate (FWER, Hochberg). If you use human ENSGs as measurements it will also map to gene names and
#' chromosome. 2) the second component returned is a measurement of "skew". It uses the absolute value of the
#' sum of the 100 largest rvals divided by the absolute value of the sum of the smallest (most negative) rvals.
#'
#'
#' @param x Matrix 1 with samples in columns and measurements
#'  (genes, transcripts, proteins, etc.), in rows. You should use
#'  the matrix from the best controlled samples for matrix 1 as
#'  it is used to identify minimally variant genes between samples.
#'  For example, if comparing differentiated neurons to human
#'  brain tissue, the differentiated neuron expression matrix
#'  should be x (matrix 1)
#' @param y Matrix 2 with same (or highly overlapping) samples
#'   and measurements as matrix 1.
#' @param method Indicate method of correlation to use "pearson", or rank order "spearman"
#' @param pct The pct of upper relative variance (variance divided by mean), to include in calculation.
#'   Many, maybe most, measurements will have minor variance that is simply noise.
#'   Measuring correlations for measurements that are relatively constant is undesirable.
#'   Set pct = .5 to only include measurements
#'   in the upper 50% of relative variance, pct = .25 to include only the upper 25%. If your data is pre-filtered
#'   to remove genes that with noise-level variance then this can be left at 1.
#' @param iter The number of iterations (permutations of x) for NULL distribution generation. If you don't want
#'   a NULL distribution analysis, set to 0
#'
#' @return
#' @importFrom dplyr %>%
#' @importFrom Hmisc rcorr
#' @export
#'
pairwiseCorSkew <- function(x, y, method = "pearson", pct = "1.0", iter = 1000)  {
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
    cat("nice matrices!\n\n")
    if (pct < 1.0) {
      cat("isolating genes in specified upper relative variance of matrix 1 measures\n\n")
      bvar <- apply(x,1,stats::var)
      bmean <- apply(x,1,base::mean)
      upv <- bvar / bmean
      combo <- cbind(x,y,upv)
      comboSort <- combo[order(-upv),]
      comboUQ <- comboSort[c(1:(nrow(comboSort)*pct)),1:(ncol(comboSort)-1)]
      x <- t(comboUQ[,c(1:(ncol(comboUQ)/2))])
      y <- t(comboUQ[,c(((ncol(comboUQ)/2)+1):ncol(comboUQ))])
    } else {
      x<- t(x)
      y <- t(y)
    }
    cat("Calculating correlations for", ncol(x), "measurements across", nrow(x), "samples\n\n")
    CWC <- NULL
    for (i in 1:ncol(x))  {
      v1 <- x[,i]
      v2 <- y[,i]
      output <- Hmisc::rcorr(v1, v2, type = method)
      Cors <- c(output$r[1,2], output$P[1,2])
      CWC <- rbind(CWC,Cors)
    }
    CWC <- cbind.data.frame(colnames(x), CWC)
    colnames(CWC) <- c("measurement","rval", "pval")
    CWC <- CWC %>%
      dplyr::mutate(FDR.BH = stats::p.adjust(pval, method = "BH", n = nrow(CWC))) %>%
      dplyr::mutate(FWER.H = stats::p.adjust(pval, method = "hochberg", n = nrow(CWC))) %>%
      dplyr::arrange(pval)

    allrhos <- CWC[,2]
    allrhos <- sort(allrhos, decreasing = TRUE)
    top100 <- allrhos[1:100]
    bottom100 <- allrhos[(length(allrhos)-99):length(allrhos)]
    skew <- abs(sum(top100)) / abs(sum(bottom100))

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
    avg <- mean(RANDrhoSkews)
    stdev <- stats::sd(RANDrhoSkews)
    zee <- (skew - avg)/stdev
    p <- stats::pnorm(q = zee, mean = 0, sd = 1, lower.tail = FALSE)

    argnames <- sys.call()
    argnames <- lapply(argnames[-1], as.character)
    arg1 <- argnames[[1]]
    arg2 <- argnames[[2]]

    grDevices::pdf(file = paste0("Correlation between ", arg1, " and ", arg2, ".pdf"))
    graphics::plot(stats::density(RANDrhoSkews), xlim = c(.4, 1.6),main = "Correlation Skew relative\n to NULL Distribution",
         xlab = "measured Skew")
    graphics::abline(v = skew, col = "red")
    graphics::abline(v = avg, col = "lightblue")
    graphics::abline(v = avg + stdev, col = "lightgrey")
    graphics::abline(v = (avg + 2*stdev), col = "lightgrey")
    graphics::abline(v = (avg + 3*stdev), col = "lightgrey")
    grDevices::dev.off()

    CWC <- list("CWC" = CWC, "skew" = skew, RANDskews = RANDrhoSkews, AvgNULL = avg, sdNULL = stdev, zScore = zee,
                pval = p, numberSamples = nrow(x), numberMeasurements = ncol(x))

    runtime <- Sys.time() - t1
    cat("It took", runtime, attr(runtime, "units"), "to run", iter, "iterations of", ncol(x), "comparisons")

    return(CWC)

  } else {print("Matrices are incompatable.\n What have you done??")}
}
