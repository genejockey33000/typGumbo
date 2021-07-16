#' @title Filter typ Object Expression Data
#'
#' @description Custom filtering of typ objects using either
#' k over a, minimum average expression, or highest x measurements
#' (based on average expression)
#'
#' @param typ.obj expression object created with 'make.typ()' function
#' @param type c("k.over.a","min.mean", ""top.x") can apply multiple
#' @param k     number of samples required to have 'a' expression level
#' @param a     expression level required
#' @param min   minimum average expression
#' @param x     number of top measurements to return
#'
#'
#' @return typ object same as entered but with filtered data and info
#'         detailing the filtering parameters
#' @export
typ.filt <- function(typ.obj, type = c("k.over.a", "min.mean", "top.x"), k = NULL, a = NULL, min = NULL, x = NULL) {
  if (!("typ" %in% class(typ.obj))) stop("This function expects an object of class'typ'")
  output <- list()
  output[["meta"]] <- typ.obj$meta
  class(output) <- class(typ.obj)
  info <- typ.obj$info
  info[["filtered"]] <- TRUE
  info[["filters"]] <- type
  filtin <- typ.obj$data

  if ("k.over.a" %in% type) {
    if (!(methods::hasArg(k))) stop("When running KoverA filtering values for k and for a must be supplied.
              k = the number of samples that must pass the expression threshold.
              a = the threshold that must be surpassed")
    if (!(methods::hasArg(a))) stop("When running KoverA filtering values for k and for a must be supplied.
              k = the number of samples that must pass the expression threshold.
              a = the threshold that must be surpassed")

      f1 <- genefilter::kOverA(k, A=a, na.rm = TRUE)
      flist <- genefilter::filterfun(f1)
      filt <- genefilter::genefilter(filtin, flist)
      filtout <- filtin[filt,]

      output[["data"]] <- filtout
      info[["KoverA.params"]] <- paste("At least",k,"samples measuring", a, "or more")
      filtin <- filtout
  }

  if ("min.mean" %in% type) {
    if (!(methods::hasArg(min))) stop("When running min.mean filtering value for 'min' must be supplied.
          min = minimum average expression across all samples.")
    chop <- apply(filtin, 1, mean) >= min
    filtout <- filtin[chop,]
    info[["Minimum.mean"]] <- min
    output[["data"]] <- filtout
    filtin <- filtout
  }

  if ("top.x" %in% type) {
    if (!(methods::hasArg(x))) stop("When running top.x filtering value for 'x' must be supplied.
                           x = number of measurements (genes or transcripts) returned ranked
                           by highest average expression")
      m <- apply(filtin, 1, mean)
      m <- cbind(m, filtin)
      m <- m[(order(m[,1], decreasing = TRUE)),]
      filtout <- m[c(1:x),-1]
      info[["filter.top.x"]] <- x
      output[["data"]] <- filtout

  }

  output[["info"]] <- info
  return(output)
}
