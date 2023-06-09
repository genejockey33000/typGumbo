#' @title Get Current Ensembl Gene Maps
#'
#' @description Returns a data frame for mapping specified genetic identifiers.
#'  For use with human data only.
#'
#' @param type Indicate type of map desired. c("full", "t2g", "t2n", "g2n", "g2e", "n2e")
#'   "full" = full map containing Ensemble transcript IDs (ENSTs), gene IDs (ENSGs), external gene names (symbols), and Entrez IDs (####)
#'     full map will contain NAs where identifiers don't have mapped counterparts
#'   "t2g" = ENST and ENSG
#'   "t2n" = ENST and gene names (symbols)
#'   "g2n" = ENSGs and gene names (symbols)
#'   "g2e" = ENSGs and Entrez IDs
#'   "n2e" = gene names (symbols) and Entrez IDs
#'
#' @return specified gene map
#' @export
#'
getEnsGeneMaps <- function(type = c("full","t2g","t2n","g2n","g2e","n2e")){
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://www.ensembl.org')
  if (type == "full") {
    map <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name", "chromosome_name"), mart = mart)
    map <- dplyr::rename(map, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, chr = chromosome_name)
    return(map)
  }
  if (type == "t2g") {
    map <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
    map <- dplyr::rename(map, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id)
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "t2n") {
    map <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
    map <- dplyr::rename(map, target_id = ensembl_transcript_id, ext_gene = external_gene_name)
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "g2n") {
    map <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name"), mart = mart)
    map <- dplyr::rename(map, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "g2e") {
    map <- biomaRt::getBM(attributes = c("ensembl_gene_id","entrezgene_id"), mart = mart)
    map <- dplyr::rename(map, ens_gene = ensembl_gene_id, entrez = entrezgene_id)
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "n2e") {
    map <- biomaRt::getBM(attributes = c("external_gene_name", "entrezgene_id"), mart = mart)
    map <- dplyr::rename(map, ext_gene = external_gene_name, entrez = entrezgene_id)
    map <- typGumbo::noNANA(map)
    return(map)
  }
}
