#' Seahorse Data import
#'
#' From mito stress test report (.xlsm) which is generated from the Wave software,
#'  copy and transpose average assay parameters calculations from the "Measures sheet"
#'  (generally starts in cell B4). Copy all measurement data up to ATP Production.
#'  Transpose data to new workbook so that sample names are in rows (merged cells) and
#'  measurements are in columns. Save as .csv file. This is the input for the seahorseImport()
#'
#' @param file quoted name of .csv file to import from above process
#'
#' @return
#' @export
seahorseImport <- function(file = NULL) {
  input <- read.csv(file = file)
  rows <- typGumbo:::odd(1:(nrow(input)))
  input2 <- input[rows,]
  row.names(input2) <- input2[,1]
  return(input2[,c(3,5,6,8)])
}
