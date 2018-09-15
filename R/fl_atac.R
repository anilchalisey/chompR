#' Create fragment size plots for ATAC-seq data
#'
#' @param bam.files character vector specifying path to pre-filtered bam files. [DEFAULT = NULL]
#' @param threads an integer value indicating the number of workers to be used. If NULL then one less than the
#' maximum number of cores will be used. [DEFAULT = NULL].
#'
#' @export
#'
#' @importFrom BiocParallel SnowParam MulticoreParam
#' @importFrom parallel detectCores
#' @importFrom csaw readParam getPESizes
#' @importFrom ggplot2 ggplot aes geom_line ggtitle theme ylab guides element_blank
#' @importFrom ggplot2 scale_x_discrete geom_hline annotate ggsave element_text xlab xlim
#' @importFrom utils write.table

fl_atac <- function(bam.files, threads = NULL) {

  if (is.null(threads)) threads <- parallel::detectCores() - 1

  if (.Platform$OS.type == "windows") {
    bparam <- BiocParallel::SnowParam(workers = threads)
  } else {
    bparam <- BiocParallel::MulticoreParam(workers = threads)
  }

  param <- csaw::readParam(pe = "both", minq = 30, dedup = TRUE, BPPARAM = bparam)

  # Get the read lengths
  run_cmd(sprintf("ls %s > samples.txt", paste(bam.files, collapse = " ")))
  rlsh <- system.file("src", "readLength.sh", package = "chompR")
  rlsh <- convert_paths(rlsh)
  rlrun_cmd <- sprintf("cat samples.txt | xargs -P %s -n 1 %s",
                       threads, rlsh)
  readlength <- as.numeric(run_cmd(rlrun_cmd, intern = T))
  unlink("samples.txt")

  # Get the fragment lengths
  PETsizes <- list()
  for (bam in bam.files) {
    PETsizes[[bam]] <-
      csaw::getPESizes(bam, param = param)
  }

  names(PETsizes) <- lapply(names(PETsizes), reduce_path)
  ps <- lapply(PETsizes, function(x) x[[1]])
  ps.plot <- list()


  for (i in seq_along(ps)) {
    ps.plot[[i]] <-
      ggplot2::ggplot() + ggplot2::aes(ps[[i]], fill = "red") +
      ggplot2::geom_density(colour = "red", alpha = 0.3) +
      ggplot2::ggtitle(paste("Fragment length distribution plot:",
                             reduce_path(names(ps)[i]))) +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1),
                     axis.ticks.x = element_blank()) +
      ggplot2::ylab("Density") + ggplot2::xlim(0, 1000) +
      ggplot2::guides(colour = FALSE, fill = FALSE) +
      ggplot2::xlab("Fragment length (bp)")
    ggplot2::ggsave(filename = paste0(
      reduce_path(bam.files)[i], "_petplot.png"),
      plot = ps.plot[[i]], device = "png")
  }

  # Write out the results
  lengths <- data.frame(row.names = reduce_path(bam.files), readlength = readlength)
  utils::write.table(lengths, "length_results.txt", row.names = TRUE,
                     col.names = TRUE, quote = FALSE, sep = "\t")
  return(lengths)
}

