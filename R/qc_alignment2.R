#' Perform alignment statistics on reads
#'
#' @param bam.files character vector specifying path to bam files
#' @param threads positive integer specifying the number of cores to use
#'
#' @export
#'
#' @importFrom ATACseqQC bamQC
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%

qc_alignment2 <- function(bam.files, threads) {

  i <- NULL

  if (threads == 1) {
    qc <- ATACseqQC::bamQC(bam.files, outPath = NULL)
    qc.sub <- qc[1:10]
    qc.df <- as.data.frame(qc.sub)
    colnames(qc.df) <- c("Reads", "Duplicated", "Mitochondrial",
                         "ProperPair", "Unmapped", "PairUnmapped",
                         "PoorQuality", "NRF", "PBC1", "PBC2")
  }

  if (threads > 1) {
    cl <- parallel::makeCluster(threads)
    doParallel::registerDoParallel(cl)
    qc <- foreach(i = 1:length(bam.files)) %dopar% {
      ATACseqQC::bamQC(bam.files[[i]], outPath = NULL)
    }
    parallel::stopCluster(cl)

    qc.sub <- lapply(qc, function(x) x[1:10])
    qc.df <- data.frame(matrix(unlist(qc.sub), nrow = length(qc.sub), byrow = TRUE))
    colnames(qc.df) <- c("Reads", "Duplicated", "Mitochondrial",
                         "ProperPair", "Unmapped", "PairUnmapped",
                         "PoorQuality", "NRF", "PBC1", "PBC2")
  }

  rownames(qc.df) <- reduce_path(bam.files)

  write.table(qc.df, file = "diagstats.txt",
              col.names = TRUE, row.names = TRUE, quote = FALSE,
              sep = "\t")

  return(qc.df)
}
