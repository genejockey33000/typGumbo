#' Gene Concept Network Plot
#' Simple generator of gene concept network plots using enrichment object
#'
#' @param x enrichment object output from running mkEnrich
#' @param nodes number of nodes (enriched pathways) to plot. default is 4
#' @param label.size How large to make the gene labels (none, teeny, small, medium, large)
#' @param category.labels LOGICAL. Should category (i.e. node or module labels) be printed (default=TRUE)
#' @param circle Logical. Should the plot be constrained to a circle (default=FALSE)
#' @param color.lines Logical. Should the lines be colored by gene set? (default=FALSE)
#'
#' @return GCN plot
#' @export
#'
gcn <- function(x, nodes=5, label.size="medium", category.labels = TRUE, circle=FALSE, color.lines=TRUE) {
    s<-.8
  if (label.size == "small") {
    s<-.5
  }
  if (label.size == "teeny"){
    s<-.3
  }
  if (label.size == "none"){
    s<-0
  }
  if (label.size == "large"){
    s<-1
  }
  if (category.labels) {
  plot <- enrichplot::cnetplot(x, showCategory = nodes,
                               cex_label_gene = s,
                               circular=circle,
                               colorEdge=color.lines)
  return(plot)} else {
    plot <- enrichplot::cnetplot(x, showCategory = nodes,
                                 cex_label_gene = s,
                                 cex_label_category = 0,
                                 circular=circle,
                                 colorEdge=color.lines)
    return(plot)
  }
}
