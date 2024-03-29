#' TSNE Plot Maker: Copied from M3C
#'
#' This is just a copy of the very nice tsne wrapper from the M3C package which uses the Rtsne
#' This is a flexible t-SNE function that can be run on a standard data frame.
#' It is a wrapper for Rtsne/ggplot2 code and can be customized with different
#' colors and font sizes and more.
#'
#' @param mydata Data frame or matrix: if dataframe/matrix should have samples as columns and rows as features
#' @param labels Character vector: if we want to just label by sex for example
#' @param perplex Numerical value: perplexity value that Rtsne uses internally
#' @param printres Logical flag: whether to print the t-SNE into current directory
#' @param seed Numerical value: optionally set the seed
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' @param textlabelsize Numerical value: text inside plot label size
#' @param legendtitle Character vector: text legend title
#' @param controlscale Logical flag: whether to control the colour scale
#' @param scale Numerical value: 1=spectral palette, 2=manual low and high palette, 3=categorical labels
#' @param low Character vector: continuous scale low color
#' @param high Character vector: continuous scale high color
#' @param colvec Character vector: a series of colors in vector for categorical labels, e.g. c("sky blue", "gold")
#' @param printheight Numerical value: png height
#' @param printwidth Numerical value: png width
#' @param text Character vector: if we wanted to label the samples with text IDs to look for outliers
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_colour_gradient
#' @importFrom ggplot2 scale_colour_distiller
#' @importFrom ggplot2  element_blank
#' @importFrom ggplot2  element_text
#' @export
tsne <- function (mydata, labels = FALSE, perplex = 15, printres = FALSE,
          seed = FALSE, axistextsize = 18, legendtextsize = 18, dotsize = 5,
          textlabelsize = 4, legendtitle = "Group", controlscale = FALSE,
          scale = 1, low = "grey", high = "red",
          colvec = c("skyblue", "gold", "violet", "darkorchid", "slateblue", "forestgreen",
                      "violetred", "orange", "midnightblue", "grey31", "black"),
          printheight = 20, printwidth = 22, text = FALSE)
{
  if (controlscale == TRUE && class(labels) %in% c("character",
                                                   "factor") && scale %in% c(1, 2)) {
    stop("when categorical labels, use scale=3")
  }
  if (controlscale == TRUE && class(labels) %in% c("numeric") &&
      scale %in% c(3)) {
    stop("when continuous labels, use scale=1 or scale=2")
  }
  if (controlscale == FALSE && scale %in% c(2, 3)) {
    warning("if your trying to control the scale, please set controlscale=TRUE")
  }
  if (sum(is.na(labels)) > 0 && class(labels) %in% c("character", "factor")) {
    warning("there is NA values in the labels vector, setting to unknown")
    labels <- as.character(labels)
    labels[is.na(labels)] <- "Unknown"
  }
  if (sum(is.na(text)) > 0 && class(text) %in% c("character", "factor")) {
    warning("there are NA values in the text vector, setting to unknown")
    text <- as.character(text)
    text[is.na(text)] <- "Unknown"
  }
  message("***t-SNE wrapper function***")
  message("running...")
  if (seed != FALSE) {
    set.seed(seed)
  }
  if (labels[1] == FALSE && text[1] == FALSE) {
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2,
                         perplexity = perplex, verbose = FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y)
    p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2)) +
      ggplot2::geom_point(colour = "skyblue", size = dotsize) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = "black"),
            axis.text.x = element_text(size = axistextsize, colour = "black"),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
      ggplot2::scale_colour_manual(values = colvec)
    if (printres == TRUE) {
      message("printing t-SNE to current directory...")
      png("TSNE.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] != FALSE && text[1] == FALSE) {
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2,
                         perplexity = perplex, verbose = FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y)
    if (controlscale == TRUE) {
      if (scale == 1) {
        p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2)) +
          ggplot2::geom_point(aes(colour = labels), size = dotsize) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize,
                colour = "black"), axis.text.x = element_text(size = axistextsize,
                colour = "black"), axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          ggplot2::labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")
      }
      else if (scale == 2) {
        p <- ggplot2::ggplot(data = scores, ggplot2::aes(x = scores$X1, y = scores$X2)) +
          ggplot2::geom_point(aes(colour =labels), size = dotsize) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                colour = "black"), axis.text.x = element_text(size = axistextsize,
                colour = "black"), axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          ggplot2::labs(colour = legendtitle) + scale_colour_gradient(low = low,
               high = high)
      }
      else if (scale == 3) {
        p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2)) +
          ggplot2::geom_point(aes(colour = labels), size = dotsize) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                colour = "black"), axis.text.x = element_text(size = axistextsize,
                colour = "black"), axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          ggplot2::labs(colour = legendtitle) + scale_colour_manual(values = colvec)
      }
    }
    else {
      p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2)) +
        ggplot2::geom_point(aes(colour = labels), size = dotsize) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize,
              colour = "black"), axis.text.x = element_text(size = axistextsize,
              colour = "black"), axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) +
        ggplot2::labs(colour = legendtitle)
    }
    if (printres == TRUE) {
      message("printing tSNE to current directory...")
      png("TSNElabeled.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] != FALSE && text[1] != FALSE) {
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2,
                         perplexity = perplex, verbose = FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y)
    scores$label <- text
    if (controlscale == TRUE) {
      if (scale == 1) {
        p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2, label = scores$label)) +
          ggplot2::geom_point(aes(colour = labels),
                    size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                    colour = "black"), axis.text.x = element_text(size = axistextsize,
                    colour = "black"), axis.title.x = element_text(size = axistextsize),
                    axis.title.y = element_text(size = axistextsize),
                    legend.title = element_text(size = legendtextsize),
                    legend.text = element_text(size = legendtextsize)) +
          ggplot2::labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral") +
          ggplot2::geom_text(vjust = "inward", hjust = "inward",
                    size = textlabelsize)
      }
      else if (scale == 2) {
        p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2, label = scores$label)) +
          ggplot2::geom_point(aes(colour = labels), size = dotsize) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize,
                colour = "black"), axis.text.x = element_text(size = axistextsize,
                colour = "black"), axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          ggplot2::labs(colour = legendtitle) +
          ggplot2::scale_colour_gradient(low = low, high = high) +
          ggplot2::geom_text(vjust = "inward", hjust = "inward", size = textlabelsize)
      }
      else if (scale == 3) {
        p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2, label = scores$label)) +
          ggplot2::geom_point(aes(colour = labels), size = dotsize) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize,
                colour = "black"), axis.text.x = element_text(size = axistextsize,
                colour = "black"), axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
          ggplot2::labs(colour = legendtitle) +
          ggplot2::scale_colour_manual(values = colvec) +
          ggplot2::geom_text(vjust = "inward", hjust = "inward", size = textlabelsize)
      }
    }
    else {
      p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2, label = scores$label)) +
        ggplot2::geom_point(aes(colour = labels), size = dotsize) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize,
              colour = "black"), axis.text.x = element_text(size = axistextsize,
              colour = "black"), axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) +
        ggplot2::labs(colour = legendtitle) +
        ggplot2::geom_text(vjust = "inward", hjust = "inward", size = textlabelsize)
    }
    if (printres == TRUE) {
      message("printing t-SNE to current directory...")
      png("TSNElabeled.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] == FALSE && text[1] != FALSE) {
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2,
                         perplexity = perplex, verbose = FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y)
    scores$label <- text
    p <- ggplot2::ggplot(data = scores, aes(x = scores$X1, y = scores$X2, label = scores$label)) +
      ggplot2::geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = dotsize) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = "black"),
            axis.text.x = element_text(size = axistextsize, colour = "black"),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
      ggplot2::scale_colour_manual(values = colvec) +
      ggplot2::geom_text(vjust = "inward", hjust = "inward", size = textlabelsize)
    if (printres == TRUE) {
      message("printing t-SNE to current directory...")
      png("TSNE.png", height = printheight, width = printwidth,
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  message("done.")
  return(p)
}
