#' Create scree plot from prcomp data
#'
#' @param pca \code{prcomp} object.
#' @inheritParams make_MA
#'
#' @importFrom ggplot2 ggplot aes geom_bar ggsave

make_scree <- function(pca) {

  scree_plot <- data.frame(pca$sdev^2/sum(pca$sdev^2))
  scree_plot[1:10,2] <- c(1:10)
  colnames(scree_plot) <- c("variance", "component number")
  scree <- ggplot2::ggplot(scree_plot[1:10,],
                           mapping = ggplot2::aes_(x = "component number", y = "variance")) +
    ggplot2::geom_bar(stat = "identity") + theme_publication()
  suppressMessages(ggplot2::ggsave('PCA_scree.png', scree, device = "png"))
}
