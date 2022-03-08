#' Clean test table
#' cleans up the sleuth_results() test table output
#'
#' @param x test table from sleuth_results() function
#' @param bcut cutoff value for fold change. Default is .5 or 50% change
#' @param qcut cutoff value for significance. Default is q value of .05
#'
#' @return
#' @export
cleanTT <- function(x, bcut = .5, qcut = .05) {
  chop <- apply(x, 1, function(x) {sum(is.na(x)) < 1})
  y <- x[chop,]
  chop <- y$ext_gene != ""
  d <- y[chop,]
  d$DE <- "NO"
  d$DE[d$b > bcut & d$qval < qcut] <- "UP"
  d$DE[d$b < -bcut & d$qval < qcut] <- "DOWN"
  d <- d %>%
    mutate(rank = -log10(qval)*b) %>%
    arrange(desc(rank))
  return(d)
}
