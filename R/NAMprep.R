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
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @export
NAMprep <- function(object) {
  allNAM.markers <- c("ROR1","GNG12-AS1","RP11-274H2.2","PLOD2","ADAMTS12","SERPINE1","CAV2","COL5A1","TEAD1","CD44","HMGA2","FRMD6","UACA","MT2A","COL6A2",
                      "PTPRC","DOCK8","SYK","APBB1IP","ALOX5AP", "ATP8B4","RUNX3","CD74","HLA-DRA",
                      "PBX1","RGS7","NRXN1","CNTNAP5","KCNH7","MAP2","DOCK3","CADM2","GAP43","ADGRL3","PPP2R2B","RIMS1","THSD7A",
                      "HECW1","PTPRN2","MIR325HG", "PAK3","CSMD1","XKR4","KCNB2","STMN2","RIMS2","CACNA1B","NCAM1","NEBL","NRG3",
                      "CNTN1","PPFIA2","FGF14","SEMA6D","RBFOX1","MAPT")
  object <- FindVariableFeatures(object)
  object <- ScaleData(object, assay = "RNA")
  scalefeats <- unique(c(object@assays$RNA@var.features, allNAM.markers))
  scalefeats <- scalefeats[scalefeats %in% row.names(object)]
  output <- Seurat::ScaleData(object, features = scalefeats)
  return(output)
}
