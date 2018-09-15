#' Generate some alignment statistics on reads
#'
#' @param samtools character string specifying path to samtools. [DEFAULT = "samtools"]
#' @param bam.files character vector specifying path to pre-filtered bam files. [DEFAULT = NULL]
#' @param filtered.bam.files character vector specifying path to filtered bam files. [DEFAULT = NULL]
#' @param threads positive integer specifying the number of cores to use
#' @param remove.mitochondrial character string.  If set, this will count reads mapping to the
#' mitochondrial genome.  The string should match the reference name for the mitochondrial genome
#' in the alignment file.  Examples include "ChrM", "M" and "MT".
#'
#' @importFrom dplyr slice
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_bar ylim theme element_rect unit labs ggtitle guides
#' @importFrom ggplot2 guide_legend geom_text coord_flip ggsave position_dodge
#'
#' @export

qc_alignment <- function(samtools = "samtools",
                         bam.files = NULL,
                         filtered.bam.files = NULL,
                         threads = 1,
                         remove.mitochondrial = "MT") {

  Sample <- value <- variable <- NULL

  outdir <- "alignmentQC"
  create_dir("alignmentQC")

  # Generate diagnostic stats from the aligned files.
  if (file.exists("temp.txt")) run_cmd("rm temp.txt")
  file.create("temp.txt", showWarnings = FALSE)

  for (i in seq_along(bam.files)) {
    code.diag <- sprintf(
      "%s view -@ %s -c %s >> temp.txt; \\
        %s view -@ %s -F 0x04 -c %s >> temp.txt; \\
        %s view -@ %s -f 1024 -c %s >> temp.txt; \\
        %s view -@ %s -c %s '%s'  >> temp.txt; \\
        %s view -@ %s -f 0x02 -c %s >> temp.txt",
      samtools, threads, bam.files[[i]], samtools, threads, bam.files[[i]],
      samtools, threads, bam.files[[i]], samtools, threads, bam.files[[i]],
      remove.mitochondrial, samtools, threads, bam.files[[i]])
    run_cmd(code.diag)
  }

  tmp <- as.vector(read.table("temp.txt"))
  diagstats <- as.data.frame(t(matrix(unlist(tmp), nrow = 5)))
  rownames(diagstats) <- basename(bam.files)
  colnames(diagstats) <- c("Total", "Mapped", "Duplicates",
                           "Mitochondrial", "ProperPair")
  diagstats$`ProperPair(%)` <- diagstats$ProperPair / diagstats$Total * 100
  diagstats$`Mapped(%)` <- diagstats$Mapped / diagstats$Total * 100
  diagstats$`Duplicates(%)` <- diagstats$Duplicates / diagstats$Total * 100
  diagstats$`Mitochondrial(%)` <- diagstats$Mitochondrial / diagstats$Total * 100
  diagstats <- diagstats[, c(1, 2, 7, 3, 8, 9, 5, 6)]

  unlink("temp.txt")

  if (!is.null(filtered.bam.files)) {
    filtered.reads <-  run_samflagstat(samtools = samtools, threads = threads,
                                       bamfile = filtered.bam.files) %>%
      dplyr::slice(1) %>% unname()
    diagstats <- as.data.frame(cbind(diagstats, t(filtered.reads)))
    colnames(diagstats)[9] <- "AfterFiltering"
    diagstats$`AfterFiltering(%)` <- diagstats$AfterFiltering / diagstats$Total * 100
  }

  # Get the read lengths
  run_cmd(sprintf("ls %s > samples.txt", paste(bam.files, collapse = " ")))
  rlsh <- system.file("src", "readLength.sh", package = "chompR")
  rlsh <- convert_paths(rlsh)
  rlrun_cmd <- sprintf("cat samples.txt | xargs -P %s -n 1 %s",
                       threads, rlsh)
  readlength <- as.numeric(run_cmd(rlrun_cmd, intern = T))
  diagstats$ReadLength <- readlength

  write.table(diagstats, file = file.path(outdir, "diagstats.txt"),
              col.names = TRUE, row.names = TRUE, quote = FALSE,
              sep = "\t")
  unlink("samples.txt")

  ## General Alignment Statistics
  if (!is.null(filtered.bam.files)) {
    df <- diagstats[, c(10, 5, 8, 3)]
    df <- data.frame("Sample" = reduce_path(rownames(df)),
                     "Remaining after filtering" = df[, "AfterFiltering(%)"],
                     "Duplicates" = df[, "Duplicates(%)"],
                     "Proper Pair" = df[, "ProperPair(%)"],
                     "Mapped" = df[, "Mapped(%)"],
                     check.names = FALSE)
  } else {
    df <- diagstats[, c(5, 8, 3)]
    df <- data.frame("Sample" = reduce_path(rownames(df)),
                     "Duplicates" = df[, "Duplicates(%)"],
                     "Proper Pair" = df[, "ProperPair(%)"],
                     "Mapped" = df[, "Mapped(%)"],
                     check.names = FALSE)
  }

  df <- reshape2::melt(df)
  df$value <- round(df$value, 1)
  p <- ggplot2::ggplot(df, ggplot2::aes(Sample, value, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", colour = "black") +
    theme_publication() + ggplot2::ylim(c(0, 110)) +
    ggplot2::theme(legend.key = ggplot2::element_rect(size = 5),
                   legend.key.size = ggplot2::unit(1.5, 'lines')) +
    ggplot2::labs(y = "Proportion (%)", x = "") +
    ggplot2::ggtitle("Alignment statistics") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = paste(value, "%")), colour = "black",
                       hjust = -0.2, position = ggplot2::position_dodge(.9), size = 3) +
    ggplot2::coord_flip()
  ggplot2::ggsave(filename = "Alignment_stats.png",
                  plot = p, device = "png",
                  path = outdir)
}
