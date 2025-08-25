#' Synaptic Vesical Release Assay Smoothing (2 Cell types)
#'
#' @param in_dir Quoted path to directory with *.xlsx file
#' @param csv LOGICAL Do you want to output *.csv files. Default is TRUE.
#' @param xlsx LOGICAL Do you want to output *.xlsx files. Default is TRUE.
#' @param filter LOGICAL Do you want to apply signal filter
#' @param split_by CHARACTER What is the name of the signal in time point 1
#' that defines the second cell type. Default is "RFP"
#' @param upper_pct_cut INTEGER upper percentage of ROIs to call positive
#' for the split_by signal. Default is 33
#'
#'@import openxlsx2
#' @returns Writes summary files for downstream usage in Prism or excel.
#' Files are organized in hierarchical directory related to sampleIDs.
#' @export
#'
syn_smooth_2color <- function(in_dir, csv = TRUE, xlsx = TRUE,
                              filter = TRUE, split_by = "RFP",
                              upper_pct_cut = 33) {
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
   Please make sure there's only one *.xlsx
   file in the directory and that it's not open. \n
   If it's closed but there's a scratch file
   (will have a ~ in front of the name), \n
   then please delete the scratch file and try again.")
  }
  if (!length(excelfile) == 1) {
    stop("Please make sure there is one and only
         one .xlsx file in the directory. Thanks!")
  }

  ## pull the full path names of all the .csv files in the in_dir
  files <- dir(in_dir, full.names = TRUE)
  excelfile <- files[grepl(pattern = ".xlsx", files)]
  wb <- openxlsx2::wb_load(excelfile)
  sheets <- openxlsx2::wb_get_sheet_names(wb)

  # trim NA's and write csv's to input directory
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
  runfiles <- paste0(sheets, ".csv")

  # make directory for output
  if (filter) {
    top_dir <- paste0(in_dir, "processed")
  } else {
    top_dir <- paste0(in_dir, "processedNoFCfilter")
  }
  dir.create(top_dir, showWarnings = FALSE)

  global_results <- wb_workbook(title = "PCT_FUN_Summary")
  ms <- NULL  #used to record number of time points for each sample

  for (f in runfiles) {
    fname <- sub("\\.csv", "", f)
    out_dirs <- c(paste0(top_dir, "/", fname, "/", split_by, "_POS"),
                  paste0(top_dir, "/", fname, "/", split_by, "_NEG"))
    full1 <- utils::read.csv(paste(in_dir, "/", f, sep = ""), header = FALSE)

    srt_row <- which(full1$V1 == "time") + 1
    srt_cols <- c(2:ncol(full1))
    srt_order <- order(as.numeric(full1[srt_row, srt_cols]), decreasing = TRUE)
    full <- full1[, -1]
    full <- full[, srt_order]
    data_cut <- which(grepl("^R.*", full[, 1]))  ## Identify row with ROI.IDs

    full_data <- full[(data_cut + 1):nrow(full), ]
    colnames(full_data) <- full[data_cut, ]
    row.names(full_data) <- NULL

    full_metrics <- full[c(1:(data_cut - 1)), ]
    row.names(full_metrics) <- full1[c(1:(data_cut - 1)), 1]
    colnames(full_metrics) <- colnames(full_data)


    # define column # cutoff
    positives <- ceiling(ncol(full_data) * (upper_pct_cut / 100))
    # subset meas to positive columns, remove upper selection column
    meas_pos <- full_data[c(2:nrow(full_data)), c(1:positives)]
    metrics_pos <- full_metrics[, c(1:positives)]
    meas_neg <- full_data[c(2:nrow(full_data)),
                          c((positives + 1):ncol(full_data))]
    metrics_neg <- full_metrics[, c((positives + 1):ncol(full_data))]

    all_meas <- list(A = meas_pos, B = meas_neg)
    all_metrics <- list(A = metrics_pos, B = metrics_neg)


    if (csv || xlsx) {
      for (out_dir in out_dirs) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      }
    }
    for (current_meas in names(all_meas)) {
      results <- wb_workbook(title = fname)
      meas <- all_meas[[current_meas]]
      meas <- apply(meas, 2, as.numeric)
      row.names(meas) <- NULL
      # define metrics matrix
      mets <- all_metrics[[current_meas]]
      mets <- t(mets)
      mets <- apply(mets, 2, as.numeric)
      row.names(mets) <- colnames(meas)

      ## Fold Change unquenched/Baseline (FC_Fun_F0) filters.  Can be
      ## toggled using filter = TRUE/FALSE argument
      fc_over_2 <- mets[, grepl("^FC_Fun.*", colnames(mets))] > 2
      fc_over_3 <- mets[, grepl("^FC_Fun.*", colnames(mets))] > 3
      total <- nrow(mets)
      pct_fc_fun_f0_greater_than_2 <- sum(fc_over_2) / total  #output
      pct_fc_fun_f0_greater_than_3 <- sum(fc_over_3) / total  #output

      # Change in absolute fluorescence stim - Baseline (DeltaF) filters
      # this filter is in place even if filter = FALSE is set
      df_over_0 <- mets[, grepl("deltaF", colnames(mets))] > 0
      pct_delta_f_greater_than_0 <- sum(df_over_0) / total  #output
      info <- rbind(pct_fc_fun_f0_greater_than_2,
                    pct_fc_fun_f0_greater_than_3,
                    pct_delta_f_greater_than_0)
      measure <- row.names(info)
      info <- cbind.data.frame(measure, info)

      if (filter) {
        cleanmeas <- meas[, fc_over_2 & df_over_0]  #output
        cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas), ]
      } else {
        cleanmeas <- meas[, df_over_0]  #output
        cleanmets <- mets[row.names(mets) %in% colnames(cleanmeas), ]
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

      # calculate mean of each time point across samples adjust mean to
      # normalized (% of max across all time points)
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
        wb_add_worksheet(sheet = "DelF") |>
        wb_add_data(x = del_f) |>
        wb_add_worksheet(sheet = "pctFun") |>
        wb_add_data(x = pct_fun) |>
        wb_add_worksheet(sheet = "info") |>
        wb_add_data(x = info)

      if (current_meas == "A") {
        active_dir <- out_dirs[grepl("POS$", out_dirs)]
        fname_sub <- paste0(fname, "_RFP_POS")
      }
      if (current_meas == "B") {
        active_dir <- out_dirs[grepl("NEG$", out_dirs)]
        fname_sub <- paste0(fname, "_RFP_NEG")
      }

      global_results <- wb_add_worksheet(global_results, sheet = fname_sub) |>
        wb_add_data(x = pct_fun)

      if (csv) {
        utils::write.csv(meas,
                         file = paste(active_dir, "/", "data.raw.csv",
                                      sep = ""), row.names = FALSE)
        utils::write.csv(mets,
                         file = paste(active_dir, "/", "metrics.raw.csv",
                                      sep = ""), row.names = TRUE)
        utils::write.csv(cleanmeas,
                         file = paste(active_dir, "/", "data.cleaned.csv",
                                      sep = ""), row.names = FALSE)
        utils::write.csv(cleanmets,
                         file = paste(active_dir, "/", "metrics.cleaned.csv",
                                      sep = ""), row.names = TRUE)
        utils::write.csv(del_f,
                         file = paste(active_dir, "/", "deltaF.csv",
                                      sep = ""), row.names = FALSE)
        utils::write.csv(pct_fun,
                         file = paste(active_dir, "/", "pctFun.csv",
                                      sep = ""), row.names = FALSE)
        utils::write.csv(info,
                         file = paste(active_dir, "/", "info.csv",
                                      sep = ""), row.names = FALSE)
      }

      if (xlsx) {
        openxlsx2::wb_workbook(title = paste0(fname_sub, "_Summary")) |>
          wb_add_worksheet(sheet = "dataRaw") |>
          wb_add_data(x = meas, row_names = FALSE) |>
          wb_add_worksheet(sheet = "metricsRaw") |>
          wb_add_data(x = mets, row_names = TRUE) |>
          wb_add_worksheet(sheet = "dataClean") |>
          wb_add_data(x = cleanmeas, row_names = FALSE) |>
          wb_add_worksheet(sheet = "metricsClean") |>
          wb_add_data(x = cleanmets, row_names = TRUE) |>
          wb_add_worksheet(sheet = "DeltaF") |>
          wb_add_data(x = del_f, row_names = FALSE) |>
          wb_add_worksheet(sheet = "pctFun") |>
          wb_add_data(x = pct_fun, row_names = FALSE) |>
          wb_add_worksheet(sheet = "info") |>
          wb_add_data(x = info, row_names = FALSE) |>
          wb_save(file = paste(active_dir, "/", fname_sub, "_Summary.xlsx",
                               sep = ""))
      }
    }
  }

  ### Following is for generating the global summary files
  ### PCT_FUN_Summary.xlsx ForPrism labelsForPrism

  # PCT_FUN_Summary
  wb_save(global_results, file = paste0(top_dir, "/PCT_FUN_Summary.xlsx"))

  ## Generate Prism friendly output
  excelfile2 <- paste0(top_dir, "/PCT_FUN_Summary.xlsx")
  pct_fun_summary <- openxlsx2::wb_load(excelfile2)
  sheets_2 <- openxlsx2::wb_get_sheet_names(pct_fun_summary)
  sheets_2_pos <- sheets_2[grepl("POS", sheets_2)]
  sheets_2_neg <- sheets_2[grepl("NEG", sheets_2)]

  ms <- max(ms) + 2

  for_prism_pos <- data.frame(matrix(ncol = 3, nrow = ms))

  for (i in sheets_2_pos) {
    temp <- openxlsx2::wb_data(pct_fun_summary, sheet = sheets_2[i])
    temp <- rbind.data.frame(colnames(temp), temp)
    temp <- temp[, c((ncol(temp) - 2):(ncol(temp)))]
    n <- rep(i, 3)
    temp <- rbind(n, temp)

    # normalize rows if unequal row numbers
    while (nrow(temp) < ms) {
      temp[nrow(temp) + 1, ] <- NA
    }

    if (i == sheets_2_pos[1]) {
      for_prism_pos[, c(1:3)] <- temp
    } else {
      for_prism_pos <- cbind.data.frame(for_prism_pos, temp)
    }
  }
  write.table(for_prism_pos,
              file = paste0(top_dir, "/RFP_POS_forPrism.csv"), sep = ",",
              row.names = FALSE, col.names = FALSE)

  ## Convert Prism labels to single character string separated by commas
  vec_row <- function(x) {
    out <- NULL

    for (i in 1:(length(x) - 1)) {
      out <- paste0(out, x[i], ",")
    }
    out <- paste0(out, x[length(x)])
    out
  }
  labels_pzm_pos <- vec_row(sheets_2_pos)
  write.table(labels_pzm_pos,
              file = paste0(top_dir, "/labelsForRFP_POS_Prism.csv"),
              quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)

  # Again for NEG
  for_prism_neg <- data.frame(matrix(ncol = 3, nrow = ms))

  for (i in sheets_2_neg) {
    temp <- openxlsx2::wb_data(pct_fun_summary, sheet = sheets_2_neg[i])
    temp <- rbind.data.frame(colnames(temp), temp)
    temp <- temp[, c((ncol(temp) - 2):(ncol(temp)))]
    n <- rep(i, 3)
    temp <- rbind(n, temp)

    # normalize rows if unequal row numbers
    while (nrow(temp) < ms) {
      temp[nrow(temp) + 1, ] <- NA
    }

    if (i == sheets_2_neg[1]) {
      for_prism_neg[, c(1:3)] <- temp
    } else {
      for_prism_neg <- cbind.data.frame(for_prism_neg, temp)
    }
  }
  write.table(for_prism_neg,
              file = paste0(top_dir, "/RFP_NEG_forPrism.csv"), sep = ",",
              row.names = FALSE, col.names = FALSE)

  ## Convert Prism labels to single character string separated by commas
  labels_pzm_neg <- vec_row(sheets_2_neg)
  write.table(labels_pzm_neg,
              file = paste0(top_dir, "/labelsForRFP_NEG_Prism.csv"),
              quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
}
