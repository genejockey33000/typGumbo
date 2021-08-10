#' Plate labels
#'
#' Builds a set of labels for diffent tissue culture plate sizes (i.e. A1, A2...). Not rocket science.
#'
#' @param plate Enter the plate type as number of wells (i.e. 6, 12, 24, 48, 96, 384)
#' @return
#' @export
platelabs <- function(plate = 96) {
  if (!(plate %in% c(6,12,24,48,96,384))) {stop("plate type should be either 6, 12, 24, 48, 96, or 384")}
  if (sum(plate %in% c(6,12,24,48,96,384)) != 1) {stop("please only enter one plate type (6, 12, 24, 48, 96, or 384)")}
  v <- NULL
  if (plate == 96) {
  rs <- c("A","B","C","D","E","F","G","H")
  cc <- 1
  c <- 12
  }
  if (plate == 48) {
    rs <- c("A","B","C","D","E","F")
    cc <- 1
    c <- 8
  }
  if (plate == 24) {
    rs <- c("A","B","C","D")
    cc <- 1
    c <- 6
  }
  if (plate == 12) {
    rs <- c("A","B","C")
    cc <- 1
    c <- 4
   }
  if (plate == 6) {
    rs <- c("A","B")
    cc <- 1
    c <- 3
  }
  if (plate == 386) {
    rs <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
    cc <- 1
    c <- 24
  }
  while (cc <= c) {
    int <- paste0(rs, cc)
    v <- c(v,int)
    cc <- cc + 1

  }

  return(v)
}
