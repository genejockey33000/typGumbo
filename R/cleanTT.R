#' Clean test table
#' cleans up the sleuth_results() test table output
#'
#' @param x test table from sleuth_results() function
#'
#' @return
#' @export
cleanTT <- function(x) {
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
