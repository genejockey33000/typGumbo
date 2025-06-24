#' TYP Proteomics Map
#' Attaches genomic mapping information to proteomics map. Generally used as an internal function in typGumbo.
#'
#' @param prot.map A simple mapping table derived from the proteomics processing that contains a uniprot ID labeled 'prot'
#' and a gene name labeled 'gene'.
#'
#' @return
#' @export
#'
#' @examples
typ.prot.map <- function(prot.map) {
  mart <- biomaRt::useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = 'https://www.ensembl.org')

  allAtts <- c("external_gene_name", "gene_biotype",
               "chromosome_name", "start_position",
               "end_position", "description")
  std.chrs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","MT")
  ens.map <- biomaRt::getBM(attributes = allAtts, mart = mart) |>
    dplyr::filter(gene_biotype == "protein_coding",
                  chromosome_name %in% std.chrs) |>
    dplyr::rename(gene = external_gene_name,
                  type = gene_biotype,
                  descr = description,
                  chr = chromosome_name,
                  start = start_position,
                  end = end_position) |>
    dplyr::select(gene, descr, chr, start, end) |>
    dplyr::mutate(prot = prot.map$prot[match(gene, prot.map$gene)], .before = 1) |>
    dplyr::filter(prot != "") |>
    dplyr::arrange(chr)
  return(ens.map)
}
