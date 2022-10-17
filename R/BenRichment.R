#' Gene Ontology Enrichment on data frame of gene symbols
#'
#'
#'
#' @param x A dataframe of gene lists in columns. Column names = name of gene list.
#'   Genes should be standardized gene symbols
#'
#' @export
BenRichment <- function(x) {
  hits <- list()
  output <- list()
  for (i in 1:(ncol(x))) {
    geneset.name <- colnames(x)[i]
    geneset <- unique(x[,i])
    geneset <- geneset[geneset != ""]
    cat("Mapping",nrow(input)," genes to entrez gene IDs\n")
    entrez <- data.frame(queryMany(geneset, scopes = 'symbol', fields = c('entrezgene'), species = 'human'), size = 1)
    hits[[geneset.name]] <- unique(entrez$entrezgene)[!(is.na(entrez$entrezgene))]
  }
  cat("\n\nRunning GO enrichment analysis for ALL ontologies... \nThis can take 1 to 2mins per geneset...\n")
  CP <- clusterProfiler::compareCluster(geneClusters = hits, fun = 'enrichGO', ont = 'ALL', OrgDb = 'org.Hs.eg.db', pAdjustMethod = 'BH', qvalueCutoff = 0.05, readable = TRUE)
  for (i in colnames(x)) {
    output[[i]] <- CP@compareClusterResult[CP@compareClusterResult$Cluster == i, ]
  }
  for (i in colnames(x)) {
    write.csv(output[[i]],file = paste0(i,".csv"), row.names = FALSE)
  }
  return(output)
}
