#' Plots of peak width and distance to TSS
#'
#' Function to plot histograms showing distribution of peak widths and distances
#' of peaks from transcriptional start sites.
#'
#' @param annotation list created by the
#' \code{annotate_peaks()} function - this contains information regarding the distance
#' from the nearest TSS.
#' @param tssdist Distance cutoff to determine proximal and distal peaks.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_histogram aes geom_line scale_x_continuous xlab
#' @importFrom ggplot2 ylab theme element_blank ggtitle annotate geom_vline aes_string
#' @importFrom GenomicRanges end start
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices png dev.off
#' @importFrom graphics hist

plot_peaksummary <- function(annotation, tssdist = 1500) {

  x <- ..density.. <- NULL

  annotated.peaks <- annotation[[1]]

  ppi <- 600
  grDevices::png(filename = paste0("peaksDistribution.png"),
                 width = 8.4*ppi, height = 6.5*ppi, res = ppi)

  widths <-
    GenomicRanges::end(annotated.peaks) - GenomicRanges::start(annotated.peaks)
  widths <- data.frame(x = widths)
  pw <-
    ggplot2::ggplot(widths, ggplot2::aes(x)) +
    ggplot2::geom_histogram(ggplot2::aes(y = (..density..)), binwidth = 100,
                            fill = "#c0392b", colour = "#c0392b", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = ..density..), stat = 'density', colour = "red") +
    ggplot2::scale_x_continuous(limits = c(0, round(max(widths), digits = -3)),
                                breaks = seq(0, round(max(widths), digits = -3),
                                             by = 500)) +
    ggplot2::xlab("peak width (bp)") + ggplot2::ylab("density") + theme_publication() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::ggtitle("Distribution of peak widths")

  tssdistances <- annotated.peaks$disttotss
  tsshist <- graphics::hist(tssdistances, plot = FALSE)
  tssdistances <- data.frame(x = tssdistances)
  tssdistances[tssdistances > 5000] <- 5000

  td <-
    ggplot2::ggplot(tssdistances, ggplot2::aes(x)) +
    ggplot2::geom_histogram(binwidth = 50, fill = "#c0392b",
                            colour = "#c0392b", alpha = 0.4) +
    theme_publication() +
    ggplot2::xlab("distance from TSS (bp)") +
    ggplot2::ggtitle("Distance of peaks from TSS") +
    ggplot2::annotate("rect", xmin = 0, xmax = tssdist, ymin = -Inf, ymax = Inf,
                      alpha = 0.1, fill = "red") +
    ggplot2::annotate("rect", xmin = tssdist, xmax = Inf, ymin = -Inf, ymax = Inf,
                      alpha = 0.1, fill = "yellow") +
    ggplot2::annotate("text", x = tssdist/2, y = max(tsshist$counts)/2,
                      label = "TSS-proximal", fontface = "bold", size = 2.2) +
    ggplot2::annotate("text", x = (tssdist + 5000)/2, y = max(tsshist$counts)/2,
                      label = "TSS-distal", fontface = "bold", size = 2.2) +
    ggplot2::geom_vline(xintercept = tssdist, lwd = 0.6, lty = 2, col = "red")

  regions <- as.data.frame(annotation[[2]][[1]])
  if (ncol(regions) == 2) colnames(regions) <- c("Region", "Frequency")
  if (ncol(regions) == 1) colnames(regions) <- c("Frequency")
  regions$Region <- c("Promoter", paste0("<", tssdist, "bp to TSS"), "5' UTR", "3' UTR", "Exon", "Intron", "Intergenic")
  regions$Region <- factor(regions$Region, levels = regions$Region)
  rp <-
    ggplot2::ggplot(regions,
                    ggplot2::aes_string(x = "Region",
                                        y = "Frequency",
                                        fill = "Region")) +
    ggplot2::geom_bar(stat = "identity", colour = "black", alpha = 0.4) +
    theme_publication() +
    ggplot2::scale_fill_manual(values = c("red", "blue", "yellow",
                                          "cyan", "orange", "pink", "grey"),
                               breaks = regions$Region,
                               labels = paste0(" ", regions$Region, "  ")) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::ylab("Percentage") +
    ggplot2::xlab("Genomic region") +
    ggplot2::ggtitle("Distribution of peaks") +
    ggplot2::theme(legend.key = ggplot2::element_rect(size = 5),
                   legend.key.size = ggplot2::unit(0.7, 'lines')) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL))

  gridExtra::grid.arrange(pw, td, rp, layout_matrix = cbind(c(1, 2), c(1, 3)))
  grDevices::dev.off()
}
