#' Add Gene Names
#'
#' Just a quick way to convert gene matrices with ENSGs as row.names to data tables with columns for human readable gene names and ENSGs.
#'
#' @param x Expression matrix (or any matrix) that has ENSG IDs as row.names
#' @param t2g Transcript to Genename table for mapping (if available). If none is supplied then will generate one from BioMart (takes longer).
#'
#' @return Data table with columns for gene names and ENSG IDs.
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @export
addGeneNames <- function(x, t2g = NULL) {
  if (is.null(t2g)) {
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://www.ensembl.org')
    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
    colnames(t2g) <- c("target_id", "ens_gene" , "ext_gene")
  }
  gene.name.map <- unique(t2g[,c("ens_gene","ext_gene")])
  ens_gene <- row.names(x)
  ext_gene <- gene.name.map$ext_gene[match(ens_gene, gene.name.map$ens_gene)]
  mat <- cbind.data.frame(ext_gene,ens_gene,x)
  return(mat)
}
