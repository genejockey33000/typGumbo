#' Synaptic Vesical Release Assay Smoothing (1 Cell type)
#'
#'
#' @param in_dir Quoted path to directory with *.xlsx file
#' @param csv LOGICAL Do you want to output *.csv files. Default is TRUE.
#' @param xlsx LOGICAL Do you want to output *.xlsx files. Default is TRUE.
#' @param filter LOGICAL Do you want to apply signal filter
#'
#' @import openxlsx2
#' @returns Writes summary files for downstream usage in Prism or excel.
#' Files are organized in hierarchical directory related to sampleIDs.
#' @export
syn_smooth_1color <- function(in_dir, csv = TRUE, xlsx = TRUE, filter = TRUE) {
  ## QC checks and cleanup in_dir entry
  if (grepl(".xlsx", in_dir) || grepl(".xls", in_dir)) {
    excelfile <- in_dir
    in_dir <- dirname(in_dir)
  } else {
    if (grepl("/$", in_dir)) {
      in_dir <- gsub("/$", "", in_dir)
    }
    files <- dir(in_dir, full.names = TRUE)
    excelfile <- files[grepl(pattern = ".xlsx", files)]
  }
  if (!(dir.exists(in_dir)))
    stop("
  \nSorry. The directory you entered doesn't seem to exist\n
                                  Please double check and try again.")

  if (sum(grepl("~", excelfile)) > 0) {
    stop("Seems like there's a scratch file open in this directory. \n
    Please make sure there's only one *.xlsx file
    in the directory and that it's not open. \n
    If it's closed but there's a scratch file
    (will have a ~ in front of the name), \n
    then please delete the scratch file and try again.")
  }

  if (!length(excelfile) == 1) {
    stop("Please make sure there is one and only one
         .xlsx file in the directory. Thanks!")
  }

  ## pull the full path names of all the .csv files in the in_dir
  files <- dir(in_dir, full.names = TRUE)
  excelfile <- files[grepl(pattern = ".xlsx", files)]
  wb <- openxlsx2::wb_load(excelfile)
  sheets <- openxlsx2::wb_get_sheet_names(wb)

  trim_mat_mess <- function(x) {
    chop_rows <- !is.na(x[, 2])
    chop_cols <- !is.na(x[1, ])
    output <- x[chop_rows, chop_cols]
    output
  }
  for (i in sheets) {
    temp <- openxlsx2::wb_read(wb, sheet = i, col_names = FALSE)
    colnames(temp) <- NULL
    temp <- trim_mat_mess(temp)
    write.csv(temp, file = paste0(in_dir, "/", i, ".csv"), row.names = FALSE)
  }

  # pull filenames for run files (csv's)
  runfiles <- paste0(sheets, ".csv")

  # make directory for output
  if (filter) {
    top_dir <- paste0(in_dir, "processed")
  } else {
    top_dir <- paste0(in_dir, "processedNoFCfilter")
  }
  dir.create(top_dir, showWarnings = FALSE)

  global_results <- wb_workbook(title = "PCT_FUN_Summary")
  ms <- NULL  #used to record number of timepoints for each sample

  for (f in runfiles) {
    fname <- sub("\\.csv", "", f)
    results <- wb_workbook(title = fname)
    out_dir <- c(paste0(top_dir, "/", fname))
    if (csv || xlsx) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    full1 <- utils::read.csv(paste(in_dir, "/", f, sep = ""), header = FALSE)
    full <- full1[, -1]
    data_cut <- which(grepl("^R.*", full[, 1]))  ## Identify row with ROI.IDs

    full_data <- full[(data_cut + 1):nrow(full), ]
    full_data <- apply(full_data, 2, as.numeric)
    colnames(full_data) <- full[data_cut, ]
    row.names(full_data) <- NULL

    full_metrics <- full[c(1:(data_cut - 1)), ]
    full_metrics <- apply(full_metrics, 2, as.numeric)
    row.names(full_metrics) <- full1[c(1:(data_cut - 1)), 1]
    colnames(full_metrics) <- colnames(full_data)

    # define measures matrix
    meas <- full_data
    # define metrics matrix
    mets <- full_metrics
    mets <- t(mets)

    ## Fold Change unquenched/Baseline (FC_Fun_F0) filters.  Can
    ## be toggled using filter = TRUE/FALSE argument
    fc_over_2 <- mets[, grepl("^FC_Fun.*", colnames(mets))] > 2
    fc_over_3 <- mets[, grepl("^FC_Fun.*", colnames(mets))] > 3
    total <- nrow(mets)
    pct_fc_fun_f0_greater_than_2 <- sum(fc_over_2) / total  #output
    pct_fc_fun_f0_greater_than_3 <- sum(fc_over_3) / total  #output

    # Change in absolute fluorescence stim - Baseline (DeltaF)
    # filters this filter is in place even if filter = FALSE is
    # set
    df_over_0 <- mets[, grepl("deltaF", colnames(mets))] > 0
    pct_delta_f_greater_than_0 <- sum(df_over_0) / total  #output
    info <- rbind(pct_fc_fun_f0_greater_than_2, pct_fc_fun_f0_greater_than_3,
                  pct_delta_f_greater_than_0)
    measure <- row.names(info)
    info <- cbind.data.frame(measure, info)

    if (filter) {
      cleanmeas <- meas[, fc_over_2 & df_over_0]  #output
      cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas),
      ]
    } else {
      cleanmeas <- meas[, df_over_0]  #output
      cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas),
      ]
    }

    # calculate measurements - average baseline (F0) for each ROI
    del_f <- NULL
    for (i in seq_len(ncol(cleanmeas))) {
      v <- cleanmeas[, i] - cleanmets[i, 1]
      del_f <- cbind(del_f, v)
    }
    colnames(del_f) <- colnames(cleanmeas)

    # calculate percentage unquenched increase over F0
    pct_fun <- NULL
    for (i in seq_len(ncol(del_f))) {
      v <- del_f[, i] / (cleanmets[i, 3] - cleanmets[i, 1])
      pct_fun <- cbind(pct_fun, v)
    }

    colnames(pct_fun) <- colnames(del_f)
    mean <- apply(pct_fun, 1, mean)
    mean_norm <- mean / max(mean)
    n <- ncol(pct_fun)
    sem <- apply(pct_fun, 1, function(x) stats::sd(x) / sqrt(n))
    pct_fun <- cbind(pct_fun, mean_norm, sem, n)
    ms <- c(ms, nrow(pct_fun))

    results <- wb_add_worksheet(results, sheet = "data.raw") |>
      wb_add_data(x = meas) |>
      wb_add_worksheet(sheet = "metrics.raw") |>
      wb_add_data(x = mets) |>
      wb_add_worksheet(sheet = "data.cleaned") |>
      wb_add_data(x = cleanmeas) |>
      wb_add_worksheet(sheet = "metrics.cleaned") |>
      wb_add_data(x = cleanmets) |>
      wb_add_worksheet(sheet = "Del_F") |>
      wb_add_data(x = del_f) |>
      wb_add_worksheet(sheet = "pctFun") |>
      wb_add_data(x = pct_fun) |>
      wb_add_worksheet(sheet = "info") |>
      wb_add_data(x = info)

    global_results <- wb_add_worksheet(global_results, sheet = fname) |>
      wb_add_data(x = pct_fun)

    if (csv) {
      utils::write.csv(meas,
                       file = paste(out_dir, "/", "data.raw.csv",
                                    sep = ""), row.names = FALSE)
      utils::write.csv(mets,
                       file = paste(out_dir, "/", "metrics.raw.csv",
                                    sep = ""), row.names = TRUE)
      utils::write.csv(cleanmeas,
                       file = paste(out_dir, "/", "data.cleaned.csv",
                                    sep = ""), row.names = FALSE)
      utils::write.csv(cleanmets,
                       file = paste(out_dir, "/", "metrics.cleaned.csv",
                                    sep = ""), row.names = TRUE)
      utils::write.csv(del_f,
                       file = paste(out_dir, "/", "deltaF.csv",
                                    sep = ""), row.names = FALSE)
      utils::write.csv(pct_fun,
                       file = paste(out_dir, "/", "pctFun.csv",
                                    sep = ""), row.names = FALSE)
      utils::write.csv(info,
                       file = paste(out_dir, "/", "info.csv",
                                    sep = ""), row.names = FALSE)
    }

    if (xlsx) {
      wb_save(results, file = paste(out_dir, "/", fname, "_Summary.xlsx",
                                    sep = ""))
    }

  }  ####END of for(f in runfiles) {} loop

  ### Following is for generating the global summary files
  ### PCT_FUN_Summary.xlsx ForPrism labelsForPrism

  # PCT_FUN_Summary
  excelfile2 <- paste0(top_dir, "/PCT_FUN_Summary.xlsx")
  wb_save(global_results, file = excelfile2)

  ## Generate Prism friendly output
  ms <- max(ms) + 2
  for_prism <- data.frame(matrix(ncol = 3, nrow = ms))

  for (i in sheets) {
    temp <- openxlsx2::wb_data(global_results, sheet = sheets[i])
    temp <- rbind.data.frame(colnames(temp), temp)
    temp <- temp[, c((ncol(temp) - 2):(ncol(temp)))]
    n <- rep(i, 3)
    temp <- rbind(n, temp)

    # normalize rows if unequal row numbers
    while (nrow(temp) < ms) {
      temp[nrow(temp) + 1, ] <- NA
    }

    if (i == sheets[1]) {
      for_prism[, c(1:3)] <- temp
    } else {
      for_prism <- cbind.data.frame(for_prism, temp)
    }
  }
  write.table(for_prism, file = paste0(top_dir, "/forPrism.csv"), sep = ",",
              row.names = FALSE, col.names = FALSE)

  ## Convert Prism labels to single character string separated by
  ## commas
  vec_row <- function(x) {
    out <- NULL
    for (i in 1:(length(x) - 1)) {
      out <- paste0(out, x[i], ",")
    }
    out <- paste0(out, x[length(x)])
    out
  }
  labels_pzm <- vec_row(sheets)
  write.table(labels_pzm, file = paste0(top_dir, "/labelsForPrism.csv"),
              quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)

}
