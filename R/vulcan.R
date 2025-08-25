#' Vulcan: Generate easy volcano plots from a "clean" test table
#'
#' Takes a clean test table and generates a volcano plot with relative ease
#' and includes a few key adjustable parameters. Test table is trimmed of all rows with
#' NAs (retained in sleuth_results()). ENSGs without mapped ext_gene names are also trimmed.
#'
#' @param x Output from CleanTable (from limma) or equivalent
#' @param labels LOGICAL. Should significant genes be labeled.
#' @param repel Force of point label repulsion (from geom_text_repel()). Default = 1
#' @param pad Padding around label boxes (from geom_text_repel()). Default = 0.25
#' @param labelsize Size of labels (from geom_text_repel()). Default = 4
#' @param title Title for volcano plot. Default is no title.
#'
#' @return volcano plot
#' @export
vulcan <- function (x, labels = TRUE, repel = 1, pad = 0.25, labelsize = 4, title) {
  d <- x
  cols <- c(UP = "#ffad73", DOWN = "#26b3ff", NO = "lightgrey")
  sizes <- c(UP = 2, DOWN = 2, NO = 1)
  alphas <- c(UP = 1, DOWN = 1, NO = 0.4)
  if (labels) {
    d$DElabel <- ""
    d$DElabel[!d$DE == "NO"] <- d$ext_gene[!d$DE == "NO"]
    p <- ggplot2::ggplot(data = d, aes(x = logFC, y = -log10(adj.P.Val),
                                       fill = DE, size = DE, alpha = DE, label = DElabel)) +
      ggplot2::geom_point(shape = 21) + ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") + ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_size_manual(values = sizes) + ggplot2::scale_alpha_manual(values = alphas) +
      ggrepel::geom_text_repel(size = labelsize, force = repel,
                               nudge_y = 0, box.padding = pad, max.overlaps = nrow(x)) +
      xlab("logFC") +
      ylab("-log10(adj.p.val)") +
      ggtitle(title)

  }
  else {
    p <- ggplot2::ggplot(data = d, aes(x = logFC, y = -log10(adj.P.Val),
                                       fill = DE, size = DE, alpha = DE)) +
      ggplot2::geom_point(shape = 21) + ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") + ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_size_manual(values = sizes) + ggplot2::scale_alpha_manual(values = alphas) +
      xlab("logFC") +
      ylab("-log10(adj.p.val)") +
      ggtitle(title)
  }
  return(p)
}
