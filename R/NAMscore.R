#' NAM (iN, Ast, iMGL) Score Generator
#'
#' Used for assigning cells to original culture when they've been mixed
#' (i.e. in triple cultures or in combined monocultures). REQUIRES that
#' ScaleData() has been applied to ALL features (not just most variable).
#' NAMprep() function can be used for that purpose if necessary.
#' NAMscore assigns an iN, Ast, and a iMGL score value for each cell.
#' Input is a Seurat object and output is a Seurat object with altered metadata
#' that includes 4 new columns. The 3 scores (iN.score, Ast.score, and iMGL.score)
#' and a "NAM.bin" column that assigns each cell to predicted culture (iN, Ast, or iMGL)
#' or as "unk". Returned Seurat object has "NAM.bin" set as active.ident
#' So running DimPlot(object) will display active reduction (probably UMAP)
#' colored by NAM.bin
#' Reset active.ident to "seurat_clusters" to display cluster info on UMAP
#' i.e.
#' object <- SetIdent(object, value = "seurat_clusters")
#'
#' Note that the markers used for Neurons, Astrocytes, and Microglia are derived from
#' snRNAseq of independent DSEB-Astros, Ngn2 induced iNs, and Blurton-Jones protocol iMGL
#' all were monocultures run on separate 10x wells thus the identities were clearly known.
#' Genetic background for all three cultures was BR24 (XX, LPNCI, APOE(E3/E3), 90yo)
#' For code generating the positive markers look in "22_0628_ipsc_runs_Dropbox/SC.SN.compare.R"
#' Strong Positives were defined be being >= 90% of cell type of interest and <= 10% of other cells.
#'
#' @param object A seurat object
#'
#' @importClassesFrom Matrix index
#' @export
#'
NAMscore <- function(object) {
  strPOS.AstIn <- c("ROR1","GNG12-AS1","RP11-274H2.2","PLOD2","ADAMTS12","SERPINE1","CAV2","COL5A1","TEAD1","CD44","HMGA2","FRMD6","UACA","MT2A","COL6A2")
  strPOS.Ast <- strPOS.AstIn[strPOS.AstIn %in% row.names(object@assays$RNA@scale.data)]
  strPOS.Ast.missing <- strPOS.AstIn[!(strPOS.AstIn %in% row.names(object@assays$RNA@scale.data))]

  strPOS.iMGLIn <- c("PTPRC","DOCK8","SYK","APBB1IP","ALOX5AP", "ATP8B4","RUNX3","CD74","HLA-DRA")
  strPOS.iMGL <- strPOS.iMGLIn[strPOS.iMGLIn %in% row.names(object@assays$RNA@scale.data)]
  strPOS.iMGL.missing <- strPOS.iMGLIn[!(strPOS.iMGLIn %in% row.names(object@assays$RNA@scale.data))]

  strPOS.iNeuroIn <- c("PBX1","RGS7","NRXN1","CNTNAP5","KCNH7","MAP2","DOCK3","CADM2","GAP43","ADGRL3","PPP2R2B","RIMS1","THSD7A",
                       "HECW1","PTPRN2","MIR325HG", "PAK3","CSMD1","XKR4","KCNB2","STMN2","RIMS2","CACNA1B","NCAM1","NEBL","NRG3",
                       "CNTN1","PPFIA2","FGF14","SEMA6D","RBFOX1","MAPT")
  strPOS.iNeuro <- strPOS.iNeuroIn[strPOS.iNeuroIn %in% row.names(object@assays$RNA@scale.data)]
  strPOS.iNeuro.missing <- strPOS.iNeuroIn[!(strPOS.iNeuroIn %in% row.names(object@assays$RNA@scale.data))]

  allNAM.markers <- c(strPOS.AstIn, strPOS.iMGLIn, strPOS.iNeuroIn)
  allNAM.markers <- allNAM.markers[allNAM.markers %in% row.names(object)]
  sub.object <- object@assays$RNA@scale.data[allNAM.markers,]

  cat(length(strPOS.Ast), " out of ",length(strPOS.AstIn), " Strong positive Ast markers detected \n")
  if (length(strPOS.Ast) != length(strPOS.AstIn)) {cat(strPOS.Ast.missing, "--not detected \n")}
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
  Ast.strength <- FeatureSetScore(sub.object = sub.object, features = strPOS.Ast)
  iMGL.strength <- FeatureSetScore(sub.object = sub.object, features = strPOS.iMGL)
  iN.strength <- FeatureSetScore(sub.object = sub.object, features = strPOS.iNeuro)

  Ast.score <- Ast.strength / (iMGL.strength + iN.strength)
  iMGL.score <- iMGL.strength / (Ast.strength + iN.strength)
  iN.score <- iN.strength / (Ast.strength + iMGL.strength)

  Ast.tru <- Ast.score > 2
  iMGL.tru <- iMGL.score > 2
  iN.tru <- iN.score > 2
  combo <- cbind(Ast.tru, iMGL.tru, iN.tru)

  guess <- rep("unk", length(Ast.score))
  guess[Ast.score > 2 & iN.score < 1 & iMGL.score < 1] <- "Ast"
  guess[iMGL.score > 2 & iN.score < 1 & Ast.score < 1] <- "iMGL"
  guess[iN.score > 2 & Ast.score < 1 & iMGL.score < 1] <- "iN"
  guess <- factor(guess, levels = c("iN", "Ast", "iMGL", "unk"))


  output[["Ast.score"]] <- Ast.score
  output[["iMGL.score"]] <- iMGL.score
  output[["iN.score"]] <- iN.score
  output[["NAM.bin"]] <- guess
  output <- SetIdent(output, value = "NAM.bin")
  return(output)
}
