#' NAM (iN, iAst, iMGL) Score Generator
#'
#' Used for assigning cells to original culture when they've been mixed
#' (i.e. in triple cultures or in combined monocultures). REQUIRES that
#' ScaleData() has been applied to ALL features (not just most variable).
#' NAMprep() function can be used for that purpose if necessary.
#' NAMscore assigns an iN, iAst, and a iMGL score value for each cell.
#' Input is a Seurat object and output is a Seurat object with altered metadata
#' that includes 4 new columns. The 3 scores (iN.score, iAst.score, and iMGL.score)
#' and a "NAM.bin" column that assigns each cell to predicted culture (iN, iAst, or iMGL)
#' or as "unk". Returned Seurat object has "NAM.bin" set as active.ident
#' So running DimPlot(object) will display active reduction (probably UMAP)
#' colored by NAM.bin
#' Reset active.ident to "seurat_clusters" to display cluster info on UMAP
#' i.e.
#' object <- SetIdent(object, value = "seurat_clusters")
#'
#' Note that the markers for iNs, iAst, and iMGLs have been hard coded and are derived from
#' scRNAseq of independent iAst, iNs, and iMGL cultures where identities were overtly known.
#' Strong Positives were defined be being >= 90% of cell type of interest and <= 10% of other cells.
#'
#' @param object A seurat object
#'
#' @return
#' @export
#'
#' @examples
  strPOS.iAstIn <- c("ROR1","GNG12-AS1","RP11-274H2.2","PLOD2","ADAMTS12","SERPINE1","CAV2","COL5A1","TEAD1","CD44","HMGA2","FRMD6","UACA","MT2A","COL6A2")
  strPOS.iAst <- strPOS.iAstIn[strPOS.iAstIn %in% row.names(object@assays$RNA@scale.data)]
  strPOS.iAst.missing <- strPOS.iAstIn[!(strPOS.iAstIn %in% row.names(object@assays$RNA@scale.data))]

  strPOS.iMGLIn <- get0(strPOS.iMGLIn, envir = asNamespace("typGumbo"))
  strPOS.iMGLIn <- c("PTPRC","DOCK8","SYK","APBB1IP","ALOX5AP", "ATP8B4","RUNX3","CD74","HLA-DRA")
  strPOS.iMGL <- strPOS.iMGLIn[strPOS.iMGLIn %in% row.names(object@assays$RNA@scale.data)]
  strPOS.iMGL.missing <- strPOS.iMGLIn[!(strPOS.iMGLIn %in% row.names(object@assays$RNA@scale.data))]

  strPOS.iNeuroIn <- get0(strPOS.iNeuroIn, envir = asNamespace("typGumbo"))
  strPOS.iNeuroIn <- c("PBX1","RGS7","NRXN1","CNTNAP5","KCNH7","MAP2","DOCK3","CADM2","GAP43","ADGRL3","PPP2R2B","RIMS1","THSD7A",
                       "HECW1","PTPRN2","MIR325HG", "PAK3","CSMD1","XKR4","KCNB2","STMN2","RIMS2","CACNA1B","NCAM1","NEBL","NRG3",
                       "CNTN1","PPFIA2","FGF14","SEMA6D","RBFOX1","MAPT")
  strPOS.iNeuro <- strPOS.iNeuroIn[strPOS.iNeuroIn %in% row.names(object@assays$RNA@scale.data)]
  strPOS.iNeuro.missing <- strPOS.iNeuroIn[!(strPOS.iNeuroIn %in% row.names(object@assays$RNA@scale.data))]

  allNAM.markers <- c(strPOS.iAstIn, strPOS.iMGLIn, strPOS.iNeuroIn)
  sub.object <- object@assays$RNA@scale.data[allNAM.markers,]

  cat(length(strPOS.iAst), " out of ",length(strPOS.iAstIn), " Strong positive iAst markers detected \n")
  if (length(strPOS.iAst) != length(strPOS.iAstIn)) {cat(strPOS.iAst.missing, "--not detected \n")}
  cat("\n",length(strPOS.iMGL), " out of ",length(strPOS.iMGLIn), " Strong positive iMGL markers detected \n")
  if (length(strPOS.iMGL) != length(strPOS.iMGLIn)) {cat(strPOS.iMGL.missing, "--not detected \n")}
  cat("\n",length(strPOS.iNeuro), " out of ",length(strPOS.iNeuroIn), " Strong positive iN markers detected \n")
  if (length(strPOS.iNeuro) != length(strPOS.iNeuroIn)) {cat(strPOS.iNeuro.missing, "--not detected \n")}

  FeatureSetScore <- function (sub.object, features) {
    m <- sub.object[features,]
    featureset.score <- 2^colSums(m)/length(features)
    return(featureset.score)
  }
  output <- object
  iAst.strength <- FeatureSetScore(sub.object = sub.object, features = strPOS.iAst)
  iMGL.strength <- FeatureSetScore(sub.object = sub.object, features = strPOS.iMGL)
  iN.strength <- FeatureSetScore(sub.object = sub.object, features = strPOS.iNeuro)

  iAst.score <- iAst.strength / (iMGL.strength + iN.strength)
  iMGL.score <- iMGL.strength / (iAst.strength + iN.strength)
  iN.score <- iN.strength / (iAst.strength + iMGL.strength)

  iAst.tru <- iAst.score > 2
  iMGL.tru <- iMGL.score > 2
  iN.tru <- iN.score > 2
  combo <- cbind(iAst.tru, iMGL.tru, iN.tru)

  guess <- rep("unk", length(iAst.score))
  guess[iAst.score > 2 & iN.score < 1 & iMGL.score < 1] <- "iAst"
  guess[iMGL.score > 2 & iN.score < 1 & iAst.score < 1] <- "iMGL"
  guess[iN.score > 2 & iAst.score < 1 & iMGL.score < 1] <- "iN"
  guess <- factor(guess, levels = c("iN", "iAst", "iMGL", "unk"))


  output[["iAst.score"]] <- iAst.score
  output[["iMGL.score"]] <- iMGL.score
  output[["iN.score"]] <- iN.score
  output[["NAM.bin"]] <- guess
  output <- SetIdent(output, value = "NAM.bin")
  return(output)
}
