#' Combine two numeric matrices with overlapping samples
#'
#' Two matrices with overlapping samples and different measurements are combined. Matrix supplied as mat1
#' is used to orient final output. Sample names must have identical formats between the two matrices and must have
#' uniquely identifiable names (i.e. "BR04" rather than vague sample names like "sample1"). That can get you into trouble
#' Samples can be in either row.names or colnames on either matrix. The script will transpose mat2 if necessary.
#' Samples not in both matrices will be removed and output will be in the same orientation (samples in rows or columns) as mat1.
#'
#' @param mat1 matrix 1 containing sample names as either row.names or colnames.
#'   Orientation of mat1 will be the orientation of the output matrix.
#' @param mat2 matrix 2 containing samples names as either row.names or colnames
#'
#' @return
#' @export
stickem <- function (mat1, mat2) {
  pare.order.rows <- function(x,y) {
    x <- x[row.names(x) %in% row.names(y),]
    y <- y[row.names(y) %in% row.names(x),]
    y <- y[match(row.names(x), row.names(y)),]
    z <- cbind(x,y)
    return(z)
  }
  pare.order.cols <- function(x,y) {
    x <- x[,colnames(x) %in% colnames(y)]
    y <- y[,colnames(y) %in% colnames(x)]
    y <- y[,match(colnames(x), colnames(y))]
    z <- rbind(x,y)
    return(z)
  }
  # ensure both matrices are oriented similarly (share rows or columns)
  a <- sum(row.names(mat1) %in% row.names(mat2))
  b <- sum(row.names(mat1) %in% colnames(mat2))
  c <- sum(colnames(mat1) %in% colnames(mat2))
  d <- sum(colnames(mat1) %in% row.names(mat2))

  if (a>b & a>c & a>d) {
    return(pare.order.rows(mat1, mat2))
  }

  if (b>a & b>c & b>d) {
    mat2 <- t(mat2)
    return(pare.order.rows(mat1, mat2))
  }

  if (c>a & c>b & c>d) {
    return(pare.order.cols(mat1,mat2))
  }

  if (d>a & d>b & d>c) {
    mat2 <- t(mat2)
    return(pare.order.cols(mat1,mat2))
  }
}
