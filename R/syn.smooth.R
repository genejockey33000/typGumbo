#' Synaptic Smoothing pipeline
#'
#' A highly specialized pipeline for processing synaptic vesicle release assay results.
#' It takes in a directory containing .csv files in a very specific TPP format (Tracy's Proprietary Processing)
#' and returns a set of processed csvs for each assay along with an excel file containing a
#' summary of the results and some analytics.
#'
#' @param in.dir directory where all input .csv files are located. Omit trailing slash.
#'     i.e. in.dir = "LViN15", or in.dir = "LViN15/results"
#' @param csv Logical. Indicate whether you want to write .csv files to output directory
#' @param xlsx Logical. Indicate whether you want to write .xlsx file with separate sheets
#'     for each output matrix to output directory. NOTE: this uses the xlsx::write.xlsx function
#'     which takes a lot of time and memory and should not be used for very large data sets (larger than
#'     1000 observations)
#' @param filter logical indicating whether data should be filtered to exclude measurements with pctFC
#'     less than 2. Even with this filter off the data will be filtered to remove measurements with pctDeltaF
#'     greater than 0 (i.e. requires positive fluorescence change)
#'
#' @export
syn.smooth <- function(in.dir, csv = TRUE, xlsx = TRUE, filter = TRUE) {
  csvs <- dir(in.dir, full.names = TRUE, recursive = TRUE)
  csvs <- csvs[grepl(pattern = ".csv", csvs)]

  runfiles <- dir(in.dir)
  runfiles <- runfiles[grepl(pattern = ".csv", runfiles)]

  if (filter) {
    out.dir <- paste0(in.dir, "processed")
  } else {
    out.dir <- paste0(in.dir, "processedNoFCfilter")
  }
  dir.create(out.dir)

  globalresults <- list()  ##added in V0.3

  for (f in runfiles) {
    results <- list()
    fname <- sub("\\.csv", "", f)
    dir.create(paste(out.dir, "/",fname, sep = ""))
    full <- utils::read.csv(paste(in.dir,"/",f, sep = ""), header = FALSE)
    full <- full[,1:(ncol(full)-4)]

    #define measurement matrix
    meas <- full[c(6:nrow(full)),c(2:ncol(full))]
    colnames(meas) <- meas[1,]
    meas <- meas[-1,]
    row.names(meas) <- c(1:nrow(meas))
    meas <- apply(meas, 2, as.numeric)


    #define metrics matrix
    mets <- full[c(1:5),]
    mets <- t(mets)
    colnames(mets) <- mets[1,]
    mets <- as.matrix(mets[-1,])
    mets <- apply(mets, 2, as.numeric)
    row.names(mets) <- colnames(meas)
    #all(colnames(meas) == row.names(mets))

    Over2 <- mets[,4] > 2
    Over3 <- mets[,4] > 3
    total <- nrow(mets)
    pct.FC_Fun_F0.greater.than.2 <- sum(Over2)/total #output
    pct.FC_Fun_F0.greater.than.3 <- sum(Over3)/total #output

    Over0 <- mets[,5] > 0
    Over1 <- mets[,5] > 1
    pct.DeltaF.greater.than.0 <- sum(Over0)/total #output
    pct.DeltaF.greater.than.1 <- sum(Over1)/total #output
    info <- rbind(pct.FC_Fun_F0.greater.than.2, pct.FC_Fun_F0.greater.than.3,pct.DeltaF.greater.than.0, pct.DeltaF.greater.than.1)
    measure <- row.names(info)
    info <- cbind.data.frame(measure, info)
    row.names(info) <- NULL

    if (filter) {
      chop <- cbind(Over2,Over0)
      chop <- apply(chop, 1, sum) == 2
      cleanmeas <- meas[,chop] #output
      cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas),]
    } else {
      cleanmeas <- meas[,Over0] #output
      cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas),]
    }

    DelF <- NULL
    for (i in c(1:(ncol(cleanmeas)))) {
      v <- cleanmeas[,i] - cleanmets[i,1]
      DelF <- cbind(DelF,v)
    }
    colnames(DelF) <- colnames(cleanmeas)

    pctFun <- NULL
    for (i in c(1:(ncol(DelF)))) {
      v <- DelF[,i]/(cleanmets[i,3]-cleanmets[i,1])
      pctFun <- cbind(pctFun,v)
    }
    colnames(pctFun) <- colnames(DelF)
    mean <- apply(pctFun, 1, mean)
    n <- ncol(pctFun)
    sem <- apply(pctFun, 1, function(x) stats::sd(x) / sqrt(n))
    pctFun <- cbind(pctFun, n, mean, sem)

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
  xlsx::write.xlsx((globalresults[[sub("\\.csv","",runfiles[1])]]), file = paste(out.dir,"/", "PCT_FUN_Summary.xlsx", sep = ""), sheetName = sub("\\.csv","",runfiles[1]), row.names = FALSE, append = FALSE)

  for (f in 2:length(runfiles)) {
    fname <- sub("\\.csv", "", f)
    xlsx::write.xlsx((globalresults[[f]]), file = paste(out.dir,"/", "PCT_FUN_Summary.xlsx", sep = ""), sheetName = sub("\\.csv","",runfiles[f]), row.names = FALSE, append = TRUE)
  }
}
