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
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @export
#'
getEnsGeneMaps <- function(type = c("full","t2g","t2n","g2n","g2e","n2e")){
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://www.ensembl.org')
  if (type == "full") {
    map <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name", "chromosome_name"), mart = mart)
    colnames(map) <- c("enst", "ensg", "gene", "chr")
        return(map)
  }
  if (type == "t2g") {
    map <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
    colnames(map) <- c("enst", "ensg")
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "t2n") {
    map <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
    colnames(map) <- c("enst", "gene")
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "g2n") {
    map <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name"), mart = mart)
    colnames(map) <- c("ensg", "gene")
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "g2e") {
    map <- biomaRt::getBM(attributes = c("ensembl_gene_id","entrezgene_id"), mart = mart)
    colnames(map) <- c("ensg", "entrez")
    map <- typGumbo::noNANA(map)
    return(map)
  }
  if (type == "n2e") {
    map <- biomaRt::getBM(attributes = c("external_gene_name", "entrezgene_id"), mart = mart)
    colnames(map) <- c("gene", "entrez")
    map <- typGumbo::noNANA(map)
    return(map)
  }
}
