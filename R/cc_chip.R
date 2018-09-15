#' Create cross-correlation plots for ChIP-seq data
#'
#' @param bam.files character vector specifying path to pre-filtered bam files. [DEFAULT = NULL]
#' @param paired logical, if TRUE then paired-end reads. [DEFAULT = TRUE]
#' @param threads an integer value indicating the number of workers to be used. If NULL then one less than the
#' maximum number of cores will be used. [DEFAULT = NULL].
#'
#' @export
#'
#' @importFrom BiocParallel SnowParam MulticoreParam
#' @importFrom parallel detectCores
#' @importFrom csaw readParam maximizeCcf correlateReads
#' @importFrom ggplot2 ggplot aes geom_line ggtitle theme ylab guides
#' @importFrom ggplot2 scale_x_discrete geom_hline annotate ggsave
#' @importFrom utils write.table

cc_chip <- function(bam.files, paired = TRUE, threads = NULL) {

  bp <- NULL

  if (is.null(threads)) threads <- parallel::detectCores() - 1

  if (.Platform$OS.type == "windows") {
    bparam <- BiocParallel::SnowParam(workers = threads)
  } else {
    bparam <- BiocParallel::MulticoreParam(workers = threads)
  }

  if (paired) {
    param <- csaw::readParam(pe = "both", minq = 30, dedup = TRUE, BPPARAM = bparam)
  } else {
    param <- csaw::readParam(pe = "none", minq = 30, dedup = TRUE, BPPARAM = bparam)
  }

  # Get the read lengths
  run_cmd(sprintf("ls %s > samples.txt", paste(bam.files, collapse = " ")))
  rlsh <- system.file("src", "readLength.sh", package = "chompR")
  rlsh <- convert_paths(rlsh)
  rlrun_cmd <- sprintf("cat samples.txt | xargs -P %s -n 1 %s",
                       threads, rlsh)
  readlength <- as.numeric(run_cmd(rlrun_cmd, intern = T))
  unlink("samples.txt")

  cc <- lapply(bam.files, function(x) csaw::correlateReads(x, param = param))
  names(cc) <- lapply(names(cc), basename)

  # turn the cross-correlation into dataframe for ggplot
  cc.df <- lapply(cc, function(x) data.frame("bp" = 1:length(x), "cc" = x))
  cc.df <- lapply(cc.df, function(x) {
    x$bp <- as.character(x$bp)
    x$bp <- factor(x$bp, levels = unique(x$bp))
    return(x)
  })


  # Get the fragment lengths
  frag.length <- lapply(cc, csaw::maximizeCcf)

  Cread <- list()
  for (i in seq_along(readlength)) {
    Cread[[i]] <- cc.df[[i]][readlength[[i]], ]$cc
  }

  Cfrag <- list()
  for (i in seq_along(cc.df)) Cfrag[[i]] <- cc.df[[i]][frag.length[[i]], ]$cc
  Cmin <- lapply(cc, min)
  cc.plot <- list()

  #calculate the NSC and RSC
  NSC <- Map("/", Cfrag, Cmin)
  RSC <- list()
  for (i in seq_along(Cfrag)) {
    RSC[[i]] <- (Cfrag[[i]] - Cmin[[i]])/(Cread[[i]] - Cmin[[i]])
  }

  for (i in seq_along(cc.df)) {
    cc.plot[[i]] <-
      ggplot2::ggplot(cc.df[[i]], ggplot2::aes(bp, cc, group = 1)) +
      ggplot2::geom_line(ggplot2::aes(colour = "red")) +
      ggplot2::ggtitle(paste("Cross-correlation plot:", reduce_path(names(cc.df)[i]))) +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)) +
      ggplot2::ylab("Cross-correlation") +
      ggplot2::guides(colour = FALSE) +
      ggplot2::scale_x_discrete(name = "Delay (bp)",
                                breaks = as.character(seq(0, 1000, 100))) +
      ggplot2::geom_hline(yintercept = Cmin[[i]],
                          colour = "blueviolet",
                          linetype = "longdash") +
      ggplot2::geom_hline(yintercept = Cread[[i]],
                          colour = "blueviolet",
                          linetype = "longdash") +
      ggplot2::geom_hline(yintercept = Cfrag[[i]],
                          colour = "blueviolet",
                          linetype = "longdash") +
      ggplot2::annotate("text", x = c(500, 500, 500),
                        y = c(Cmin[[i]] + Cfrag[[i]]*0.02,
                        Cread[[i]] + Cfrag[[i]]*0.02,
                        Cfrag[[i]]*1.02),
                        label = c("Cmin", "Cread", "Cfrag"),
                        size = 3) +
      ggplot2::annotate("text", x = 180,
                        y = Cmin[[i]] + Cfrag[[i]]*0.02,
                        fontface = "bold",
                         label = paste0("NSC = ", round(NSC[[i]], 2),
                                        "    RSC = ", round(RSC[[i]], 2)))
    ggplot2::ggsave(filename = paste0(
      reduce_path(bam.files)[i], "_ccplot.png"),
      plot = cc.plot[[i]], device = "png")
  }

  # Write out the results
  lengths <- data.frame(row.names = reduce_path(bam.files), readlength = readlength,
                        fragmentlength = unlist(frag.length), NSC = unlist(NSC), RSC = unlist(RSC))
  utils::write.table(lengths, "length_results.txt", row.names = TRUE,
                     col.names = TRUE, quote = FALSE, sep = "\t")
  return(lengths)
}

