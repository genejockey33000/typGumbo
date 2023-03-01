#' Synaptic Smoothing pipeline
#'
#' A highly specialized pipeline for processing synaptic vesicle release assay results.
#' It takes in a directory containing an excel file (.xlsx) in a very specific TPP format (Tracy's Proprietary Processing)
#' and returns a set of processed csvs for each assay along with an excel file containing a
#' summary of the results and some analytics. Each sheet from the input excel file should have 5 header rows from row 1 they
#' must be 1) F0 (average starting fluorescence), 2) Fstim (average stim fluorescence), 3) Fun (max unquenched fluro), 4) FC_Fun-F0 (fold change of
#' unquenched max relative to baseline(F0)), and 5) deltaF (increase in fluor in stimulated vs. F0). First column (under metrics labels)
#' must be a space (extra column) for time values. The script will remove this column and if it contains data rather than
#' time stamps then the data will be lost. The worksheets used to generate the csv's will often have extra columns at the end of the sheet.
#' They start at row 6 (after the header metrics) and report average and error fluor across all measurements for a given time. There can be
#' 2-4 of these or none. It doesn't matter as any column after the last value in the row 1 metrics will be removed. If this function
#' returns a Java out of memory heap space error then try restarting R and then run the following line before loading typGumbo
#' options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
#'
#'
#'
#' @param in.dir directory where all input .csv files are located.
#'     i.e. in.dir = "LViN15", or in.dir = "LViN15/results". Can also enter path directly
#'     to excel file i.e. in.dir = "LViN15/LViN15.xlsx"
#' @param csv Logical. Indicate whether you want to write .csv files to output directory
#' @param xlsx Logical. Indicate whether you want to write .xlsx file with separate sheets
#'     for each output matrix to output directory. NOTE: this uses the xlsx::write.xlsx function
#'     which takes a lot of time and memory and should not be used for very large data sets (larger than
#'     1000 observations)
#' @param filter logical. indicating whether data should be filtered to exclude measurements with pctFC
#'     less than 2. Even with this filter off the data will be filtered to remove measurements with pctDeltaF
#'     greater than 0 (i.e. requires positive fluorescence change)
#'
#' @importFrom readxl excel_sheets
#' @export
syn.smooth <- function(in.dir, csv = TRUE, xlsx = TRUE, filter = TRUE) {
  ## pull the full path names of all the .csv files in the in.dir
  if (grepl(".xlsx", in.dir) | grepl(".xls", in.dir)) {
    excelfile <- in.dir
    in.dir <- dirname(in.dir)
  } else {
    if (grepl("/$", in.dir)) {
      in.dir <- gsub("/$", "", in.dir)}
    files <- dir(in.dir, full.names = TRUE)
    excelfile <- files[grepl(pattern = ".xlsx", files)]
  }
  if (!(dir.exists(in.dir))) stop("
  \nso ummm.  here's the thing... \nthat directory doesn't actually exist on your machine. \nmaybe you put a leading slash in front of a relative path?\n
I mean I don't know what you did but I don't think that directory is there. \n
Try again?\n
With an actual existing path?")
  sheets <- readxl::excel_sheets(excelfile)

  trimMatMess <- function(x) {
    chopRows <- !is.na(x[,2])
    chopCols <- !is.na(x[1,])
    return(x[chopRows[,1],chopCols[1,]])
  }
  for ( i in sheets ) {
    temp <- readxl::read_xlsx(path = excelfile, sheet = i, col_names = FALSE)
    colnames(temp) <- NULL
    temp <- trimMatMess(temp)
    write.csv(temp, file = paste0(in.dir,"/",i, ".csv"), row.names = FALSE)
  }
  csvs <- dir(in.dir, full.names = TRUE, recursive = TRUE)
  csvs <- csvs[grepl(pattern = ".csv", csvs)]

  # pull filenames for run files (csv's)
  runfiles <- dir(in.dir)
  runfiles <- runfiles[grepl(pattern = ".csv", runfiles)]

  # make directory for output
  if (filter) {
    out.dir <- paste0(in.dir, "processed")
  } else {
    out.dir <- paste0(in.dir, "processedNoFCfilter")
  }
  dir.create(out.dir)

  globalresults <- list()  ##added in V0.3
  Ms <- NULL

  for (f in runfiles) {
    results <- list()
    fname <- sub("\\.csv", "", f)
    #read and trim data from csv
    dir.create(paste(out.dir, "/",fname, sep = "")) #create sub-directory for individual csv
    full <- utils::read.csv(paste(in.dir,"/",f, sep = ""), header = FALSE) #read current csv

    #define measurement matrix
    meas <- full[c(6:nrow(full)),c(2:ncol(full))] #remove metrics rows (top 5) and time stamps (first column)
    colnames(meas) <- meas[1,]
    meas <- meas[-1,]
    row.names(meas) <- c(1:nrow(meas))
    meas <- apply(meas, 2, as.numeric) #convert from character to numeric


    #define metrics matrix
    mets <- full[c(1:5),]
    mets <- t(mets)
    colnames(mets) <- mets[1,]
    mets <- as.matrix(mets[-1,])
    mets <- apply(mets, 2, as.numeric)
    row.names(mets) <- colnames(meas)
    #all(colnames(meas) == row.names(mets))  #should return true (for testing)

    ## Fold Change unquenched/Baseline (FC_Fun_F0) filters. Can be toggled using filter = TRUE/FALSE argument
    FC_Over2 <- mets[,4] > 2
    FC_Over3 <- mets[,4] > 3  ##currently unused
    total <- nrow(mets)
    pct.FC_Fun_F0.greater.than.2 <- sum(FC_Over2)/total #output
    pct.FC_Fun_F0.greater.than.3 <- sum(FC_Over3)/total #output

    # Change in absolute fluorescence stim - Baseline (DeltaF) filters
    # this filter is in place even if filter = FALSE is set
    DF_Over0 <- mets[,5] > 0
    pct.DeltaF.greater.than.0 <- sum(DF_Over0)/total #output
    info <- rbind(pct.FC_Fun_F0.greater.than.2, pct.FC_Fun_F0.greater.than.3,pct.DeltaF.greater.than.0)
    measure <- row.names(info)
    info <- cbind.data.frame(measure, info)
    row.names(info) <- NULL

    if (filter) {
      chop <- cbind(FC_Over2,DF_Over0)
      chop <- apply(chop, 1, sum) == 2
      cleanmeas <- meas[,chop] #output
      cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas),]
    } else {
      cleanmeas <- meas[,DF_Over0] #output
      cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas),]
    }

    # calculate measurements - average baseline (F0) for point of interest
    DelF <- NULL
    for (i in c(1:(ncol(cleanmeas)))) {
      v <- cleanmeas[,i] - cleanmets[i,1]
      DelF <- cbind(DelF,v)
    }
    colnames(DelF) <- colnames(cleanmeas)

    # calculate percentage unquenched increase over F0
    pctFun <- NULL
    for (i in c(1:(ncol(DelF)))) {
      v <- DelF[,i]/(cleanmets[i,3]-cleanmets[i,1])
      pctFun <- cbind(pctFun,v)
    }
    colnames(pctFun) <- colnames(DelF)
    mean <- apply(pctFun, 1, mean)
    n <- ncol(pctFun)
    sem <- apply(pctFun, 1, function(x) stats::sd(x) / sqrt(n))
    pctFun <- cbind(pctFun, mean, sem, n)
    Ms <- c(Ms, nrow(pctFun))

    results[["data.raw"]] <- meas
    results[["metrics.raw"]] <- mets
    results[["data.cleaned"]] <- cleanmeas
    results[["metrics.cleaned"]] <- cleanmets
    results[["DelF"]] <- DelF
    results[["pctFun"]] <- pctFun
    results[["info"]] <- info

    globalresults[[fname]] <- pctFun

    if (csv) {
      utils::write.csv(meas, file = paste(out.dir,"/",fname,"/","data.raw.csv", sep = ""), row.names = FALSE)
      utils::write.csv(mets, file = paste(out.dir,"/",fname,"/","metrics.raw.csv", sep = ""), row.names = TRUE)
      utils::write.csv(cleanmeas, file = paste(out.dir,"/",fname,"/","data.cleaned.csv", sep = ""), row.names = FALSE)
      utils::write.csv(cleanmets, file = paste(out.dir,"/",fname,"/","metrics.cleaned.csv", sep = ""), row.names = TRUE)
      utils::write.csv(DelF, file = paste(out.dir,"/",fname,"/","deltaF.csv", sep = ""), row.names = FALSE)
      utils::write.csv(pctFun, file = paste(out.dir,"/",fname,"/","pctFun.csv", sep = ""), row.names = FALSE)
      utils::write.csv(info, file = paste(out.dir,"/",fname,"/","info.csv", sep = ""), row.names = FALSE)
    }

    if (xlsx) {
      xlsx::write.xlsx(meas, file = paste(out.dir,"/",fname,"/",fname,"summary.xlsx", sep = ""), sheetName = "dataRaw", row.names = FALSE, append = FALSE)
      xlsx::write.xlsx(mets, file = paste(out.dir,"/",fname,"/",fname,"summary.xlsx", sep = ""), sheetName = "metricsRaw", row.names = TRUE, append = TRUE)
      xlsx::write.xlsx(cleanmeas, file = paste( out.dir,"/",fname,"/",fname,"summary.xlsx", sep = ""), sheetName = "dataClean", row.names = FALSE, append = TRUE)
      xlsx::write.xlsx(cleanmets, file = paste(out.dir,"/",fname,"/",fname,"summary.xlsx", sep = ""), sheetName = "metricsClean", row.names = TRUE, append = TRUE)
      xlsx::write.xlsx(DelF, file = paste(out.dir,"/",fname,"/",fname,"summary.xlsx", sep = ""), sheetName = "DeltaF", row.names = FALSE, append = TRUE)
      xlsx::write.xlsx(pctFun, file = paste(out.dir,"/",fname,"/",fname,"summary.xlsx", sep = ""), sheetName = "pctFun", row.names = FALSE, append = TRUE)
      xlsx::write.xlsx(info, file = paste(out.dir,"/",fname,"/",fname,"summary.xlsx", sep = ""), sheetName = "info", row.names = FALSE, append = TRUE)
    }

  }
  fname <- sub("\\.csv", "", runfiles[1])
  xlsx::write.xlsx((globalresults[[fname]]), file = paste(out.dir,"/", "PCT_FUN_Summary.xlsx", sep = ""), sheetName = fname, row.names = FALSE, append = FALSE)

  for (f in 2:length(runfiles)) {
    fname <- sub("\\.csv", "", runfiles[f])
    xlsx::write.xlsx((globalresults[[fname]]), file = paste(out.dir,"/", "PCT_FUN_Summary.xlsx", sep = ""), sheetName = fname, row.names = FALSE, append = TRUE)
  }
  ## Generate Prism friendly output
  Ms <- max(Ms) + 2
  forPrism <- data.frame(matrix(ncol = 3, nrow = Ms))
  excelfile2 <- paste0(out.dir, "/PCT_FUN_Summary.xlsx")
  sheets2 <- readxl::excel_sheets(path = excelfile2)
  for ( i in sheets2 ) {
    temp <- readxl::read_xlsx(path = paste(out.dir,"/", "PCT_FUN_Summary.xlsx", sep = ""), sheet = i, col_names = FALSE)

    temp <- temp[,c((ncol(temp)-2):(ncol(temp)))]

    n <- rep(i, 3)
    temp <- rbind(n, temp)

    #normalize rows if unequal row numbers
    while (nrow(temp) < Ms) {
      temp[nrow(temp) + 1,] <- NA
    }

    if (i == sheets2[1]) {
      forPrism[,c(1:3)] <- temp
    } else {
      forPrism <- cbind.data.frame(forPrism, temp)
    }
  }
  write.table(forPrism, file = paste0(out.dir, "/forPrism.csv"), sep = ",",row.names = FALSE, col.names = FALSE)

  ## Convert Prism labels to single character string separated by commas
  VecRow <- function(x) {
    out <- NULL

    for(i in 1:(length(x)-1)) {
      out <- paste0(out, x[i], ",")
    }
    out <- paste0(out,x[length(x)])
    return(out)
  }

  labelsPzm <- VecRow(sheets2)
  write.table(labelsPzm, file = paste0(out.dir, "/labelsForPrism.csv"), quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
}
