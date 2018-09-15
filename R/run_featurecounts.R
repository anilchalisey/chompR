#' Count reads in bam files using featureCounts
#'
#' @description  Function to count reads mapping to user-provided regions.
#'
#' @param featurecounts character string specifying path to featureCounts. [DEFAULT = "featurecounts"]
#' @param annotationFile character string, specifying path to region file in SAF format.
#' @param requireBothEndsMapped logical, if TRUE, only fragments that have both
#'   ends successfully aligned will be considered for summarization. This option
#'   should be used together with pairedEnd = TRUE. [DEFAULT = TRUE]
#' @param excludeChimeric logical, if TRUE, the chimeric fragments (those
#'   fragments that have their two ends aligned to different chromosomes) will
#'   NOT be counted. This option should be used together with pairedEnd = TRUE. [DEFAULT = TRUE]
#' @param pairedEnd logical, if TRUE,  fragments (or templates) will be counted
#'   instead of reads. This option is only applicable for paired-end reads. [DEFAULT = TRUE]
#' @param fragmentLength vector of two integers specifying the minimum and
#'   maximum fragment/template lengths.  Must be used together with pairedEnd = TRUE.
#'   [DEFAULT = c(50, 600)].
#' @param countMultiMapping logical, if TRUE then multi-mapping reads/fragments will be
#'   counted. [DEFAULT = FALSE]
#' @param multiFeatureReads logical, if TRUE then reads/fragments overlapping with more than one
#'   meta-feature/feature will be counted more than once. Note that when
#'   performing meta-feature level summarization, a read (or fragment) will
#'   still be counted once if it overlaps with multiple features within the same
#'   meta-feature (as long as it does not overlap with other metafeatures). [DEFAULT = FALSE]
#' @param minQual integer, the minimum mapping quality score a read must satisfy in
#'   order to be counted. For paired-end reads, at least one end should satisfy
#'   this criteria. [DEFAULT = 0]
#' @param stranded Indicate if strand-specific read counting should be
#'   performed.  Acceptable values: 0 (unstranded), 1 (stranded) and 2
#'   (reversely stranded). For paired-end reads, strand of the
#'   first read is taken as the strand of the whole fragment. [DEFAULT = 0]
#' @param threads an integer value indicating the number of workers to be used. If NULL then one less than the
#' maximum number of cores will be used. [DEFAULT = NULL].
#' @param ignoreDup logical, if TRUE then reads that were marked as duplicates will
#'   be ignored.  In paired end data, the entire read pair will be ignored if at
#'   least one end is found to be a duplicate read. [DEFAULT = FALSE]
#' @param outname Character string.  Name of the output file. The output file
#' contains the number of reads assigned to each meta-feature or feature.
#' @param alignments character vector specifying path to alignment files. [DEFAULT = NULL]
#' @param quiet logical, if \code{TRUE} then summary information not displayed
#'
#' @return
#' Tab-delimited file containing the number of reads assigned to each
#' meta-feature or feature and a matrix with the same.
#'
#' @importFrom data.table fwrite fread
#' @importFrom parallel detectCores
#' @importFrom utils head tail
#'
#' @export
#'
#' @examples
#' \dontrun{
#' alignments <- list.files(path = "bam_files", pattern = "*.bam$")
#' threads <- parallel::detectCores() - 1
#' counts <- count_reads(featurecounts = "featureCounts", threads = threads,
#'                       alignments = alignments)
#' }

run_featurecounts <- function(featurecounts = "featureCounts",
                              annotationFile = NULL,
                              requireBothEndsMapped = TRUE,
                              excludeChimeric = TRUE,
                              pairedEnd = TRUE,
                              fragmentLength = c(30, 800),
                              countMultiMapping = FALSE,
                              multiFeatureReads = FALSE,
                              minQual = 0,
                              stranded = 0,
                              threads = NULL,
                              ignoreDup = FALSE,
                              outname = NULL,
                              alignments = NULL,
                              quiet = TRUE) {


  if (is.null(outname)) outname <- "raw_counts.txt"
  if (is.null(threads)) threads <- parallel::detectCores() - 1

  feature.counts <- paste(featurecounts,
                          "-a", annotationFile,
                          "-F SAF",
                          if (isTRUE(requireBothEndsMapped)) {"-B"},
                          if (isTRUE(excludeChimeric)) {"-C"},
                          if (isTRUE(pairedEnd)) {"-p"},
                          if (!is.null(fragmentLength)) {
                            paste("-P","-d", fragmentLength[1],
                                  "-D", fragmentLength[2])},
                          if (isTRUE(countMultiMapping)) {"-M"},
                          if (isTRUE(multiFeatureReads)) {"-O"},
                          "-Q", minQual,
                          "-s", stranded,
                          "-T", threads,
                          if (isTRUE(ignoreDup)) {"--ignoreDup"},
                          "-o", paste0(outname, ".full"),
                          paste(alignments, collapse = " "))
  give_note(feature.counts)
  res <- run_cmd(feature.counts)
  results <- as.data.frame(data.table::fread(input = paste0(outname, ".full"),
                                             skip = 1))
  results <- results[, -c(2:6)]
  colnames(results) <- reduce_path(colnames(results))
  rownames(results) <- results$Geneid
  results <- as.matrix(results[, -1])

  if (!quiet) {
    give_note(paste("\nTotal number of features:\n", nrow(results)))
    Sys.sleep(1)
    give_note("\nHead of counts matrix:\n")
    print(utils::head(results))

    Sys.sleep(1)
    give_note("\nTail of counts matrix:\n")
    print(utils::tail(results))
  }

  write.table(results, file = outname, col.names = NA, row.names = TRUE,
              sep = "\t", quote = FALSE)

  return(results)
}
