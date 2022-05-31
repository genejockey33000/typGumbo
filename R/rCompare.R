#' Overlay Waterfall plot of pos and negative correlations
#' Takes output of pairwiseCorrSkew() and generates an inverted
#' waterfall plot visually comparing pos correlations, FDR passers
#' and neg correlations. Plot is outputted as a pdf
#'
#' @param x Output list object from running pairwiseCorrSkew
#' @param title Title of plot (describe the comparison). I.e. "iAstro to Brain"
#' @param pdfName Quoted name of the pdf file to generate as output
#'
#' @return
#' @export
#'
rCompare <- function(x, title="add title", pdfName="output.pdf") {
  mat <- x$CWC
  mcps <- mat[mat$FDR.BH <= 0.05,]
  pos <- mat[mat$rval >= 0,]
  neg <- mat[mat$rval < 0,]
  lin <- function(x) {
    y <- x$rval
    names(y) <- x$measurement
    return(y)
  }

  mcps <- lin(mcps)
  pos <- lin(pos)
  neg <- lin(neg)
  neg <- sort(neg, decreasing = TRUE)
  negFlp <- sort(abs(neg), decreasing = TRUE)

  xxis1 <- sum(length(pos), length(neg))
  xxis2 <- max(c(length(pos), length(neg)))

  maxr <- max(c(pos, negFlp))
  yxis <- maxr*1.1

  no.mpcs <- length(mcps)

  pdf(file = pdfName)
  par(mar = c(5, 5, 5, 2))
  plot(pos,
                  pch=".",
                  type = "h",
                  col = "aquamarine3",
                  main = paste0("All Correlations \n", title),
                  ylim = c(-yxis,yxis),
                  ylab = "R value",
                  xlim = c(0, xxis1),
                  xlab = "Genes",
                  cex.lab = 1.7,
                  cex.main = 1.8
  )
  lines(pos, col = "grey50")
  lines(c(rep(0, (length(pos))),neg), col = "grey50")
  lines(c(rep(0,(length(pos))),neg), type = "h", col = "lightsalmon")
  lines(pos[1:(no.mpcs)], type = "h", col = "turquoise4")
  legend(x=(xxis1 * .6),y=maxr, fill = c("turquoise4", "aquamarine3", "lightsalmon"), legend = c( "Pos. FDR <.05", "Pos. Correlations", "Neg. Correlations"))
  abline(h=0, col = "darkgrey")

  par(mar = c(5, 5, 5, 2))
  plot(pos,  #plot positives and set axes
                 pch=".",
                 type = "h",
                 col = "aquamarine3",
                 main = paste0("All Correlations \n", title),
                 ylim = c(0,.9),
                 ylab = "abs(R value)",
                 xlab = "Genes",
                 cex.lab = 1.7,
                 cex.main = 1.8
  )
  lines(pos, col = "grey50") #outline positives
  lines(c(rep(0,(no.mpcs + 1)),negFlp), col = "grey50") #outline negs starting after end of FDR passers
  lines(c(rep(0,(no.mpcs + 1)),negFlp), type = "h", col = "lightsalmon")  #fill in negatives
  lines(pos[1:(no.mpcs)], type = "h", col = "turquoise4") #fill in FDR passers
  legend(x=(xxis2 * .6),y=maxr, fill = c("turquoise4", "aquamarine3", "lightsalmon"), legend = c( "Pos. FDR <.05", "Pos. Correlations", "Neg. Correlations"))
  dev.off()
}
