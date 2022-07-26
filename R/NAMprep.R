#' NAM prep
#' Running the NAMscore function on a Seurat object uses the scale.data slot
#' If Scale() was only run on the most variable features (default) rather than
#' all features then many of the key markers could be lost (especially if one
#' cell type is over-represented). NAM prep process a Seurat object to run the
#' ScaleData() function on all features rather than only the 2000 most variable.
#' Can take a while.
#'
#' @param object Seurat object
#'
#' @return
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @export
NAMprep <- function(object) {
  load("R/sysdata.rda")
  scalefeats <- unique(c(object@assays$RNA@var.features, allNAM.markers))
  scalefeats <- scalefeats[scalefeats %in% row.names(object)]
  output <- Seurat::NormalizeData(object) %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData(features = scalefeats) %>%
    Seurat::RunPCA() %>%
    Seurat::RunUMAP(dims = 1:30) %>%
    Seurat::FindNeighbors(k.param = 10, dims = 1:30) %>%
    Seurat::FindClusters()
  return(output)
}
