#' Run sleuth on kallisto output
#'
#' For somewhat easier processing of kallisto output. If you have your kallisto output files
#' in a sub-directory named 'results' and a matching metadata file named 'meta' then you only
#' need to execute 'output <- run.sleuth() to output both gene-level and transcript-level
#' sleuth results (returned as a list). Can specify  custom directory or metadata file names,
#' and specify "gene" or "transcript" levels in arguments.
#'
#'
#' @param in.dir Name of directory with kallisto output. Sub direstories should be structured as:
#' sampleName1/{abundance.h5,abundance.tsv,run_info.json}
#' sampleName2/{abundance.h5,abundance.tsv,run_info.json}
#' ...
#'
#' @param meta Name of meta data file in working directory i.e. "AH1-18meta.csv"
#' @param level Level of analysis to perform "gene" or "trans". Default is both.
#'
#' @return sleuth object(s)
#' @export
run.sleuth <- function(in.dir = "results", meta = "meta.csv", level = c("gene")) {
  if (level != "gene" & level != "trans") stop('Please set level = "gene" or level = "trans".')
  path <- file.path(in.dir, dir(in.dir))
  s2c <- utils::read.csv(file = paste(meta), header = TRUE, stringsAsFactors = FALSE)
  colnames(s2c)[1] <- "sample"
  s2c <- cbind.data.frame(s2c,path)

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
  colnames(t2g) <- c("target_id", "ens_gene" , "ext_gene")

  if ("gene" %in% level) {
  sogene <- sleuth::sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE, target_mapping = t2g, gene_mode = TRUE, aggregation_column = 'ens_gene')
  return(sogene)
  }
  if ("trans" %in% level) {
  sotrans <- sleuth::sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE, target_mapping = t2g)
  return(sotrans)
  }
}
