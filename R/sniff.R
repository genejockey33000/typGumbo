#' @title Quick peek (sniff) of typ objects
#'
#' @description  Just a quick peek. Nothing complicated
#'
#' @param typ.obj object created with typGumbo::make.typ() function
#' @param what name of field to display. Default is what = "meta" options are "gene" (snippet
#'     of gene-level matris), "trans" (snippet of transcription-level matrix), "info"
#'     (processing information), "all" (all of the above).
#'
#' @return Quick View
#' @export
sniff <- function(typ.obj, what = "meta") {
  if (!("typ" %in% class(typ.obj))) stop("This function expects an object of class'typ'")
  if ("gene" %in% what | "all" %in% what) {
    gene <- typ.obj$datagene
  utils::View(gene[1:20,])
  }
  if ("trans" %in% what | "all" %in% what) {
    trans <- typ.obj$datatrans
  utils::View(trans[1:20,])
  }
  if ("meta" %in% what | "all" %in% what) {
    meta <- typ.obj$meta
  utils::View(meta)
  }
  if ("info" %in% what | "all" %in% what) {
    info <- typ.obj$info
  utils::View(info)
  }
}
