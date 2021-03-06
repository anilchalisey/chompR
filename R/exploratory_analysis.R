#' Exploratory analysis of DESeq2 data
#'
#' @param dds DESeqDataSet object created by \code{DESeq2::DESeq}
#' @inheritParams deseq2_analysis
#'
#' @importFrom grDevices png dev.off
#' @importFrom DESeq2 plotDispEsts vst
#' @importFrom SummarizedExperiment assay
#' @importFrom utils write.table

exploratory_analysis <- function(dds, species, metadata) {
  grDevices::png(
    filename = "DESeq2_dispersionplot.png",
    width = 8,
    height = 8,
    units = "in",
    res = 500)
  DESeq2::plotDispEsts(dds, genecol = "grey",
                       fitcol = "purple",
                       finalcol = "orange")
  grDevices::dev.off()

  vst <- DESeq2::vst(dds)
  normcounts <- SummarizedExperiment::assay(vst)
  rlddf <- data.frame(normcounts)
  rownames(rlddf) <- names(dds)
  colnames(rlddf) <- dds$sample

  make_PCA(normcounts, metadata$sampleinfo)

  utils::write.table(
    x = rlddf,
    file = "DESeq2_vst_normalisedcounts.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}
