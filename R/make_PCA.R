#' Create PCA plot for transformed data
#'
#' @param normcounts normalised count data; created as part of
#' \code{edger_analysis}, \code{limma_voom_analysis}, or \code{deseq2_analysis}.
#' @param sample.info metadata created by \code{sanity_check}
#' @inheritParams make_MA
#'
#' @importFrom matrixStats rowVars
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab ggsave
#' @importFrom ggrepel geom_text_repel

make_PCA <- function(normcounts, sample.info) {

  PC1 <- PC2 <- Condition <- NULL

  rv <- matrixStats::rowVars(normcounts)
  pca <- stats::prcomp(t(normcounts[order(rv, decreasing = TRUE)[
    seq_len(min(500, length(rv)))
    ], ]))
  make_scree(pca)
  dpt <- data.frame(Condition = sample.info$condition,
                    sample = sample.info$sample, pca$x[, 1:2])
  varprop <- summary(pca)$sdev^2
  varprop <- varprop/sum(varprop)
  varprop <- paste0(round(varprop[1:2]*100, 1), "%")

  pcaplot <- ggplot2::ggplot(data = dpt, ggplot2::aes(x = PC1, y = PC2, label = sample)) +
    ggplot2::geom_point(ggplot2::aes(colour = Condition), size = 3) +
    ggrepel::geom_text_repel() +
    ggplot2::xlab(paste("PC1:", varprop[1], "of variance")) +
    ggplot2::ylab(paste("PC2:", varprop[2], "of variance")) +
    theme_publication()
  pname.file <- "PCAplot_normalizedcounts.png"
  suppressMessages(ggplot2::ggsave(pname.file, pcaplot))
}


