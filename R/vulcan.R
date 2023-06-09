#' Vulcan: Generate easy volcano plots from sleuth results test tables
#'
#' Takes the output of sleuth_results() and generates a volcano plot with relative ease
#' and includes a few key adjustable parameters. Test table is trimmed of all rows with
#' NAs (retained in sleuth_results()). ENSGs without mapped ext_gene names are also trimmed.
#'
#' @param x Output from sleuth_results() or any test tables with values for b and for qval
#' @param bcut Cutoff value for b (log fold change) for significance
#' @param qcut Cutoff value for qval (adjusted pval) for significance
#' @param labels How many of the top DEG data points to label. Top DEGs are identified
#'   buy calculating the rank = -lot10(qval) * bval. Vulcan will plots the top and bottom scores.
#'   Default is 0 and prints no labels.
#' @param repel Force of point label repulsion (from geom_text_repel()). Default = 1
#' @param pad Padding around label boxes (from geom_text_repel()). Default = 0.25
#' @param labelsize Size of labels (from geom_text_repel()). Default = 4
#'
#' @return volcano plot
#' @export
vulcan <- function(x, bcut = .5, qcut = .05, labels = 0, repel = 1, pad = .25, labelsize = 4) {
  d <- cleanTT(x, bcut = bcut, qcut = qcut)
  cols <- c("UP" = "#ffad73", "DOWN" = "#26b3ff", "NO" = "lightgrey")
  sizes <- c("UP" = 2, "DOWN" = 2, "NO" = 1)
  alphas <- c("UP" = 1, "DOWN" = 1, "NO" = 0.4)

  if (labels > 0) {
    d$DElabel <- NA
    d$DElabel[1:labels] <- d$ext_gene[1:labels]
    d$DElabel[(nrow(d)-(labels-1)):nrow(d)] <- d$ext_gene[(nrow(d)-(labels-1)):nrow(d)]

    p <- ggplot2::ggplot(data=d, aes(x=b, y=-log10(qval),
                            fill = DE,
                            size = DE,
                            alpha = DE,
                            label = DElabel)) +
      ggplot2::geom_point(shape = 21) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_size_manual(values = sizes) +
      ggplot2::scale_alpha_manual(values = alphas) +
      ggrepel::geom_text_repel(size = labelsize, force = repel, nudge_y = 5, box.padding = pad)
  } else {
    p <- ggplot2::ggplot(data=d, aes(x=b, y=-log10(qval),
                            fill = DE,
                            size = DE,
                            alpha = DE)) +
      ggplot2::geom_point(shape = 21) +
      ggplot2::theme_minimal() +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_size_manual(values = sizes) +
      ggplot2::scale_alpha_manual(values = alphas)
  }
  return(p)
}
