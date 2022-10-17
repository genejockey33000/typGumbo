#' All Duplicated
#' Return logical vector marking any value that has a duplication of another value.
#'  Not just the duplicates of smaller indexes.
#'
#' @param vec Input vector that might have duplicate values
#'
#' @export
allDuplicated <- function(vec){
  front <- duplicated(vec)
  back <- duplicated(vec, fromLast = TRUE)
  all_dup <- front + back > 0
  return(all_dup)
}
