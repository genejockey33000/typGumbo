#' Quick PCA
#'
#' Just an abbreviated wrapper for prcomp() and fviz_pca_ind because I get tired of all the lines of code. Cleans up scripts with many PCA plots.
#'
#' @param x Your expression matrix with samples in columns and measurements (genes/proteins/etc.) in rows.
#' @param meta Metadata object, ideally with column labeled "sample" that matches the names of the columns in your expression matrix. all(colnames(x) == meta[,sample]) must return TRUE
#' @param sample Name of metadata column with sample names. Default is "sample" specify if yours is different (i.e. "Sample", "samp." or some other foolish name)
#' @param scale (Logical) Whether or not prcomp should scale the data prior to running PCA. Leave FALSE if data is already scaled/log transformed. Set to TRUE for raw data values
#' @param colorBy (Optional) Name of the column to color dots by. (i.e. "treatment", "batch", "sex")
#' @param samp.labels (Logical) Should individual points be
#' @param plot.title (Optional) Title for your PCA plot.
#'
#' @return PCA plot
#' @importFrom stats prcomp
#' @importFrom factoextra fviz_pca_ind
#'
#' @export
quickPCA <- function(x, meta = NULL, sample = "sample", scale = FALSE, colorBy = NULL, samp.labels = FALSE, plot.title = NULL) {
  if(all(colnames(x) == meta[,sample])) {

    noVar <- apply(x, 1, function(x){var(x) == 0}) # test for 0 variance measurements and remove if necessary
    if (sum(noVar) > 0) {x <- x[!noVar,]}

    input.pca <- stats::prcomp(base::t(x), scale. = scale)
    if (samp.labels) {
      plot <- factoextra::fviz_pca_ind(input.pca,
                   axes = c(1,2),
                   col.ind = as.factor(meta[,colorBy]),
                   geom = c("point", "text"),
                   pointsize = 3,
                   palette = "rbPal",
                   addEllipses = FALSE,
                   legend.title = colorBy,
                   repel = TRUE
      )} else {
        plot <- factoextra::fviz_pca_ind(input.pca,
                     axes = c(1,2),
                     col.ind = as.factor(meta[,colorBy]),
                     geom = c("point"),
                     pointsize = 3,
                     palette = "rbPal",
                     addEllipses = FALSE,
                     legend.title = colorBy,
                     repel = TRUE
        )
      }
    return()
  } else{stop("The samples (column names) in the matrix don't match the samples in the metadata")}
}
