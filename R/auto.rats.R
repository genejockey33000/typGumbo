#' @title Autosomal Ratios (auto.rats)
#'
#' @description For each sample, calculates ratio of normalized
#'    expression for all measurements from each autosome to all
#'    other autosomes (i.e. autosome ratio). Then compares ratios
#'    across all samples for each chromosome. Can detect severe
#'    chromosomal abnormalities (i.e. chromosomal duplication or loss).
#'    Outputs data frame of all ratios with column for each chromosome
#'    and a pdf containing all ratio plots for chr1:chr22. Pdfs are
#'    written to a sub-folder named 'plots'
#'
#' @param typ.obj typ objects of type 'typ.gene' or 'typ.trans'
#' @param color.by OPTIONAL Specify column in object meta data
#'                 to color bar plots (i.e. "sex" or "karyotype")
#'
#' @return matrix of ratios, pdf containing plots
#' @export
auto.rats <- function(typ.obj, color.by = NULL) {
  if (!("typ" %in% class(typ.obj))) stop("This function expects an object of class'typ'")
  if ("prot" %in% class(typ.obj)) stop("I don't think this should be run on protein data.
            BUT, if YOU think there's a reason it should, let me know and maybe I can tweak this
            script to allow that. ;-) -Richard")

  df <- typ.obj$data
    if ("gene" %in% class(typ.obj)) {
      mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
      mapper <- biomaRt::getBM(attributes = c("ensembl_gene_id","chromosome_name"), mart = mart)
      colnames(mapper) <- c("ens_gene", "chr")
    ens_gene <- row.names(df)
    df <- cbind.data.frame(ens_gene,df)
    row.names(df) <- NULL
    df <- merge(mapper, df, by.x = "ens_gene", by.y = "ens_gene", all.x = FALSE, all.y = TRUE)
    chop <- is.na(df$chr)
    df <- df[!chop,]
    } else {
      mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
      mapper <- biomaRt::getBM(attributes = c("ensembl_transcript_id","chromosome_name"), mart = mart)
      colnames(mapper) <- c("ens_trans", "chr")

      ens_trans <- row.names(df)
      df <- cbind.data.frame(ens_trans,df)
      row.names(df) <- NULL
      df <- merge(mapper, df, by.x = "ens_trans", by.y = "ens_trans", all.x = FALSE, all.y = TRUE)
      chop <- is.na(df$chr)
      df <- df[!chop,]
    }

  autosomes <- as.character(1:22)
  df <- df[(df$chr %in% autosomes),]

  if (methods::hasArg(color.by)) {
  output <- typ.obj$meta[,c("sample", color.by)]
  } else { output <- typ.obj$meta[,"sample"] }

  for (cc in autosomes) {
    n.cc <- autosomes[(autosomes) != cc]
    cc.df <- df[df$chr == cc,]
    n.cc.df <- df[df$chr %in% n.cc,]
    avg.cc <- apply(cc.df[,(colnames(df) %in% typ.obj$meta[,"sample"])],2,mean)
    avg.ncc <- apply(n.cc.df[,(colnames(df) %in% typ.obj$meta[,"sample"])],2,mean)
    ratio <- avg.cc/avg.ncc
    output <- cbind.data.frame(output, ratio)
    colnames(output)[(ncol(output))] <- paste("chr",cc, sep = "")
    colnames(output)[1] <- "sample"
  }

  chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
  plots <- list()
  maindir <- getwd()
  plotdir <- "plots"
  dir.create(file.path(maindir, plotdir), showWarnings = FALSE)

  if (methods::hasArg(color.by)) {
  grDevices::pdf("plots/AutosomeRatioPlots.pdf")
  for (chr in chrs) {
    thisPlot <- ggplot2::ggplot(output, ggplot2::aes(x = stats::reorder(sample, !!ggplot2::ensym(chr)), fill=!!ggplot2::ensym(color.by), y = !!ggplot2::ensym(chr))) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_brewer(palette = "Dark2") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 90)) +
      ggplot2::ggtitle(paste(chr,"rel to other autosomes")) +
      ggplot2::xlab("ratio") +
      ggplot2::ylab("sample")

    plots[[chr]] <- thisPlot
    print(thisPlot)
  }
  grDevices::dev.off()
  }  else {
  grDevices::pdf("plots/AutosomeRatioPlots.pdf")
  for (chr in chrs) {
    thisPlot <- ggplot2::ggplot(output, ggplot2::aes(x = stats::reorder(sample, !!ggplot2::ensym(chr)), y = !!ggplot2::ensym(chr))) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_brewer(palette = "Dark2") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 90)) +
      ggplot2::ggtitle(paste(chr,"rel to other autosomes")) +
      ggplot2::xlab("ratio") +
      ggplot2::ylab("sample")

    plots[[chr]] <- thisPlot
    print(thisPlot)
  }
  grDevices::dev.off()

  }
  return(output)
}
