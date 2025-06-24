#' Plot Cluster Composition
#'
#' Takes a Seurat object as input and uses embedded metadata to generate a set of stacked bar plots illustrating percentages of
#'
#' @param so a Seurat object
#' @param cluster
#' @param comp.report
#' @param table.out
#'
#' @return
#' @export
#'
#' @examples
plot_cluster_composition = function (so, cluster = NULL, comp.report = NULL, table.out = FALSE) {
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)

  count_table <- table(unlist(so@meta.data[cluster]), unlist(so@meta.data[comp.report]))
  count_mtx   <- as.data.frame.matrix(count_table)
  total <- rowSums(count_mtx)
  out <- cbind(count_mtx, total)
  if(table.out) {return(out)}

  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  sorted_labels <- paste(sort(levels(cluster_size$cluster),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  colnames(melt_mtx)[2] <- "dataset"


  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") +
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
    scale_fill_brewer(palette = "Spectral") +
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")

  p2 + p1 + plot_layout(widths = c(2.5,1.5))
}
