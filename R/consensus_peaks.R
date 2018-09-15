#' Consensus peaks
#'
#' Given a GRangesList with peaks for each sample, this function determines
#' consensus peaks (found in at least N replicates, where N is defined by the
#' user)
#'
#' @param peaks A GRangesList object consisting of one GRanges object for each
#'              set of peaks.
#' @param replicates The minimum number of replicates in which a peak must be
#'                   present for it to count as a consensus peak.
#' @param plot Boolean to indicate whether or not to generate diagnostic plots.
#'             Default = FALSE.
#'
#' @return
#' A GRanges object of consensus peaks
#'
#' @importFrom ggplot2 ggplot geom_bar aes guides theme element_blank ylab ggtitle geom_text
#' @importFrom GenomicRanges reduce findOverlaps elementMetadata
#' @importFrom grDevices png dev.off
#'
#' @export
#'
#' @examples
#' \dontrun{
#' conPeaks <- consensus_peaks(peaks = peaks, replicates = 2)
#' }

consensus_peaks <- function(peaks, replicates, plot = FALSE) {
  # Concatenate all peaks and then reduce them so overlapping peaks are merged
  allpeaks <- peaks[[1]]
  for (i in seq_along(peaks)[-1]) {
    allpeaks <- c(allpeaks, peaks[[i]])
  }
  allpeaks <- GenomicRanges::reduce(allpeaks)

  # For each reduced peak, determine whether it was present in each sample type
  pr <- lapply(peaks, function(x) GenomicRanges::findOverlaps(allpeaks, x))
  dr <- data.frame(matrix(nrow = length(allpeaks), ncol = length(peaks),
                          dimnames = list(1:length(allpeaks))))
  colnames(dr) <- names(peaks)
  for (i in seq_along(dr)) {
    dr[, i][as.numeric(rownames(dr)) %in% pr[[i]]@from] <- 1
    dr[, i][!(as.numeric(rownames(dr)) %in% pr[[i]]@from)] <- 0
  }
  dr$All <- apply(dr, MARGIN = 1, function(x) sum(x))
  GenomicRanges::elementMetadata(allpeaks) <- dr

  # Find regions that are present in at least 'n' replicates
  consensuspeaks <- allpeaks[allpeaks$All >= replicates]

  if (isTRUE(plot)) {

    # Some basic stats
    df <- as.data.frame(c(sapply(peaks, length), allpeaks = length(allpeaks),
                          consensus = length(consensuspeaks)))
    colnames(df) <- "No of peaks"
    print(df)

    df$sample <- as.character(rownames(df))
    df$sample <- factor(df$sample, levels = unique(df$sample))

    ppi <- 600
    grDevices::png(filename = "consensus.png",
                   width = 8.4*ppi,
                   height = 6.5*ppi,
                   res = ppi)
    bp <- ggplot2::ggplot(df, ggplot2::aes(x = df$sample, y = df$`No of peaks`)) +
      ggplot2::geom_bar(stat = "identity",
                        ggplot2::aes(fill = df$sample), colour = "black") +
      theme_publication() + ggplot2::guides(fill = FALSE) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
      ggplot2::ylab("No. of peaks") +
      ggplot2::ggtitle("Peak counts per sample and for merged and consensus tracks") +
      ggplot2::geom_text(ggplot2::aes(label = df$`No of peaks`),
                         vjust = -0.5, fontface = "bold")
    print(bp)
    grDevices::dev.off()
  }
  return(consensuspeaks)
}
