#' Run DESeq2 analysis on ChIP/ATAC-seq count data
#'
#' @inheritParams run_chip
#' @param metadata list created by \code{sanity_check_chip} and with raw
#' count data stored in \code{metadata$counts}.  Created as part of
#' \code{run_chip}
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom DESeq2 DESeqDataSetFromMatrix counts DESeq results
#' @importFrom edgeR cpm
#' @importFrom utils combn write.table
#' @importFrom grDevices png dev.off
#' @importFrom dplyr as_tibble arrange desc
#'
#' @export

deseq2_analysis <- function(metadata, species) {

  log2FoldChange <- padj <- NULL

  threads <- metadata$threads
  doParallel::registerDoParallel(threads)

  deseqdata <- DESeq2::DESeqDataSetFromMatrix(
    countData = metadata$counts,
    colData = metadata$sampleinfo,
    design = metadata$design
  )

  minLibSize <- min(colSums(DESeq2::counts(deseqdata)))
  minGroupSize <- min(tabulate(deseqdata$condition))
  dds <- deseqdata
  dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = TRUE)

  create_dir(file.path(metadata$outdir, "DESeq2"))
  out <- file.path(metadata$outdir, "DESeq2")
  tryCatch(exploratory_analysis(dds, species = species, metadata = metadata))
  tomove <- list.files(".", pattern = "DESeq2_")
  move_file(tomove, out)

  contrasts <- factor(levels(dds$condition), levels = levels(metadata$sampleinfo$condition))
  con <- utils::combn(contrasts, 2)
  vec <- list()
  for (i in 1:ncol(con)) {
    vec[[i]] <- levels(con[, i, drop = TRUE])
  }

  contrasts <- lapply(vec, function(x) c("condition", rev(x)))
  res.names <- lapply(contrasts, function(x) paste(x[2], x[3], sep = "-"))
  out.sub <- file.path(out, res.names)
  lapply(out.sub, create_dir)

  DEG <- list()
  for (i in seq_along(contrasts)) {
    res <- DESeq2::results(dds, parallel = TRUE,
                           contrast = contrasts[[i]], independentFiltering = FALSE)

    columns.of.interest = c("peak", "baseMean", "log2FoldChange", "pvalue", "padj")
    output <- cbind(peak = row.names(res), as.data.frame(res), stringsAsFactors = FALSE)[, columns.of.interest]

    make_MA(output, fdr = 0.05, label.rectangle = TRUE,
            top = 20, select.method = "logfc")
    move_file(tomove = "DESeq2_MAplot.png", out.sub[[i]])

    DEG[[i]] <- subset(output, padj < 0.05)
    utils::write.table(
      x = as.data.frame(output),
      file = file.path(out.sub[[i]], "DESeq2_differential_expression.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE)
    }

  names(DEG) <- res.names

  lapply(names(DEG), function(x) {
    return(
      give_note(paste("\nFound", nrow(subset(DEG[[x]], log2FoldChange > 0)), "upregulated peak(s) and",
                      nrow(subset(DEG[[x]], log2FoldChange < 0)), "downregulated peaks(s) for the contrast", x,
                      "using DESeq2.\n\n", collapse=" "))
    )
  })

  dev.off()

  return(lapply(DEG, function (x) dplyr::as_tibble(x) %>%
                dplyr::arrange(dplyr::desc(abs(log2FoldChange)), padj)))
}
