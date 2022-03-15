#' Get Current Ensembl Gene Maps
#'
#' Requires no input. Returns a list object with 2 data frames for mapping ENSGs and ENSTs to gene symbols.
#'  For use with human genomic data only.
#'
#'  tgn data frame has columns for
#'  transcript ID (ENST numbers) labeled "target_id"
#'  gene ID (ENSG numbers) labeled "ens_gene"
#'  gene symbols (human readable gene symbols) labeled "ext_gene
#'
#'  gn contains only ens_gene, and ext_gene columns (unique)
#'
#' @return
#' @export
#'
getEnsGeneMaps <- function(){
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
  tgn <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)
  tgn <- dplyr::rename(tgn, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  gn <- unique(tgn[,c(2,3)])
  output <- list(tgn = tgn, gn = gn)
  return(output)
}
