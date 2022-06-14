#' Cluster-wise test for DEGs of Seurat Object
#'
#' For analysis within each cluster (or other set of cell group defined in meta.data)
#' for DEGs (markers) that are significant between specified comparison.
#' Comparison must exist as a column in seuratObject(at)meta.data" column with
#' two possible values (i.e. ko, wt, or treated, sham)
#'
#' @param x Seurat object with clusters defined
#' @param group Category to compare with ("seurat_clusters" is default). Can be
#'   any parameter specified in meta.data column that segregated cells into labeled groups.
#' @param comp Comparison to be tested. Must be specified and be a column in meta.data.
#'   Use quoted column header of comparison, i.e. comp = "treatment" where "treatment" is a
#'   column header for cells that were "stim." or "vehicle".
#'
#' @return
#'
#' @importFrom Seurat FindMarkers
#' @export
#'
DEbyClust <- function(x, group = "seurat_clusters", comp = NULL) {
  if(group != "seurat_clusters"){
    clusters <- unique(x@meta.data[,group])
  } else {clusters <- levels(x@meta.data$seurat_clusters)}

  c1 <- unique(x@meta.data[,comp])[1]
  c2 <- unique(x@meta.data[,comp])[2]
  output <- list()
  DE.sum <- NULL
  DE.all <- NULL
  n <- 1

  for(c in clusters) {
    oldw <- getOption("warn")
    options(warn = -1)

    id1 <- row.names(x@meta.data)[x@meta.data[,group] == c & x@meta.data[,comp] == c1]
    id2 <- row.names(x@meta.data)[x@meta.data[,group] == c & x@meta.data[,comp] == c2]
    DE <- Seurat::FindMarkers(x, ident.1 = id1, ident.2 = id2)

    if(nrow(DE) != 0){

      name <- paste0("Clust.", c, ".", comp, ".", c1, ".v.", c2)
      DE <- cbind.data.frame(c, comp, DE)
      feature <- row.names(DE)
      DE <- cbind.data.frame(feature, DE)
      row.names(DE) <- NULL
      colnames(DE)[1:3] <- c("feature", group, "compare")

      DE.sig <- DE[(DE[,"p_val_adj"] < 0.05),]

      if(nrow(DE.sig) == 0) {
        DE.sig <- DE[1,]
      }

      DE.sum <- rbind.data.frame(DE.sum, DE.sig)
      DE.all <- rbind.data.frame(DE.all, DE)
      output[[name]] <- DE
    }
    cat("completed", n, "comparisons", length(clusters)-n, "to go")
    n<-n+1
  }
  options(warn = oldw)
  output[["sig.summary"]] <- DE.sum
  output[["all.summary"]] <- DE.all
  return(output)
}
