#' Plot Cluster Composition
#'
#' Takes a Seurat object as input and uses embedded metadata to generate a
#' set of stacked bar plots illustrating compositional percentages
#'
#' @param so a Seurat object
#' @param cluster Cluster number or ID to analyse
#' @param comp.report name of the metadata parameter to report composition
#' @param table.out LOGICAL. Should only table values be reported
#'
#' @export
#'
plot_cluster_composition = function (so, cluster = NULL, comp.report = NULL, table.out = FALSE) {
  library(patchwork)
  library(ggplot2)
  library(reshape2)

  count_table <- table(unlist(so@meta.data[cluster]), unlist(so@meta.data[comp.report]))
  count_mtx   <- as.data.frame.matrix(count_table)
  total <- rowSums(count_mtx)
  out <- cbind(count_mtx, total)
  if(table.out) {return(out)}

  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- reshape2::melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- stats::aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  sorted_labels <- paste(sort(levels(cluster_size$cluster),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  colnames(melt_mtx)[2] <- "dataset"


  p1 <- ggplot2::ggplot(cluster_size, aes(y= cluster,x = value)) +
    ggplot2::geom_bar(position="dodge", stat="identity",fill = "grey60") +
    ggplot2::theme_bw() + ggplot2::scale_x_log10() +
    ggplot2::xlab("Cells per cluster, log10 scale") + ggplot2::ylab("")
  p2 <- ggplot2::ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
    ggplot2::geom_bar(position="fill", stat="identity") +
    ggplot2::theme_bw() + ggplot2::coord_flip() +
    ggplot2::scale_fill_brewer(palette = "Spectral") +
    ggplot2::ylab("Fraction of cells in each dataset") +
    ggplot2::xlab("Cluster number") +
    ggplot2::theme(legend.position="top")

  p2 + p1 + patchwork::plot_layout(widths = c(2.5,1.5))
}
