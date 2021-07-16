#' @title Make Standardized 'typ' object
#'
#' @description Generates uniform R object of class 'typ' for use in
#'    downstream processes
#'
#' @param x sleuth object or expression matrix
#' @param  level    c("gene", "trans", "prot")
#' @param  meta    meta data information as data frame (not required if sleuth obj passed)
#' @param  regress OPTIONAL column names in meta data file specifying regression
#'
#' @return object of class 'typ'
#' @export
make.typ <- function(x, level = "gene", meta = NULL, regress = NULL) {
  output <- list()
  if (!("sleuth" %in% class(x))) {
  class(output) <- c("typ",level)
  }
  info <- list()
  if (length(regress) > 3) stop("make.typ only allows 3 regression parameters")
  if (length(level) > 1) stop("Each typ object should only have one dataset ('gene' or 'trans' or 'prot').
                              You can combine typ.objects into 'collections' downstream.")
  if (!("sleuth" %in% class(x)) && !methods::hasArg(meta)) stop("If you don't supply a sleuth object then you must supply a meta data file (as a dataframe please)")
  if (!("sleuth" %in% class(x)) && !all(regress %in% colnames(meta))) stop("One or more of your batches aren't in your metadata")
  if (!("sleuth" %in% class(x))) {info[["source"]] <- "supplied.matrix"}

  if ("sleuth" %in% class(x)) {
    info[["source"]] <- "sleuth.object"
    meta <- as.data.frame(x$sample_to_covariates)
    param <- colnames(meta)
    if (methods::hasArg(regress) && !(all(regress %in% colnames(meta)))) stop("At least one of your regression parameters doesn't match a metadata column name.
    Check available options using
    colnames({sleuth.object.name}$sample_to_covariates)
    for sleuth objects
    or
    colnames(meta)
    for manual metadata")
    if (x$gene_mode == "TRUE") {level = "gene"} else {level = "trans"}
    class(output) <- c("typ",level)
    x <- sleuth::sleuth_to_matrix(obj = x, which_df = "obs_norm", which_units = "tpm")
    info[["source.df"]] <- "obs_norm"
    info[["source.units"]] <- "tpm"
  }
  info[["level"]] <- level

    if (!("sample" %in% colnames(meta))) stop("metadata file must have column labeled 'sample'")
    if (!(all(meta$sample == colnames(x)))) stop("sample column in meta needs to match column names in x")
    if (level == "gene" | level == "trans") {
      lowv <- apply(x, 1, stats::var) != 0
      x <- x[lowv, ]
      if (level == "gene") {ens_gene <- verOut(row.names(x))
      row.names(x) <- ens_gene
      } else {ens_trans <- verOut(row.names(x))
      row.names(x) <- ens_trans}
      x <- as.matrix(x)

      if (methods::hasArg(regress)) {
        CBin <- log2(x + 0.1)
        cbmod <- stats::model.matrix(~1, data = meta)
        batch1 <- as.factor(meta[, regress[1]])
        CBout <- sva::ComBat(dat = CBin, batch = batch1, mod = cbmod)
        if (length(regress) > 1) {
          batch2 <- as.factor(meta[, regress[2]])
          CBout <- sva::ComBat(dat = CBout, batch = batch2, mod = cbmod)
          if (length(regress) > 2) {
            batch3 <- as.factor(meta[, regress[3]])
            CBout <- sva::ComBat(dat = CBout, batch = batch3, mod = cbmod)
          }
        }
        output[["data"]] <- 2^(CBout)

      } else {output[["data"]] <- x}

    } else {
      if (level != "prot") stop("level must be set to either 'gene', 'trans', or 'prot'")
      lowv <- apply(x, 1, stats::var) != 0
      x <- x[lowv, ]
      x <- as.matrix(x)
      if (methods::hasArg(regress)) stop("I'm sorry but I can't regress protein data during the typ object build just yet.
                                Pester Richard about this if you need it.")
      output[["data"]] <- x
    }

  output[["meta"]] <- meta
  output[["info"]] <- info

  return(output)
}
