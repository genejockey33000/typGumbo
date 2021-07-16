#' Sleuth object to Matrix (internal function)
#'
#' Exports sleuth object to an expression matrix using normalized expression
#'     values and expressed as transcripts per million (TPM) units
#'
#' @param obj
#' @param which_df
#' @param which_units
#'
function (obj, which_df, which_units) {
  data <- as.data.frame(obj[[which_df]])
  res <- list()
  s_data <- dplyr::select_(data, "target_id", "sample", which_units)
  s_data <- tidyr::spread_(s_data, "sample", which_units)
  rownames(s_data) <- s_data$target_id
  s_data$target_id <- NULL
  s_data <- as.matrix(s_data)
  s_data
}
