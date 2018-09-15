#' Import the metadata and check it is correct
#'
#' @param sample.info character string giving the path to a tab-delimited text
#' file with at least the columns <condition> (treatment condition),
#' <sample> (sample name), and <file1> (absolute or relative path to the
#' fastq files).  If fastq files and PE reads, then a column
#' <file2> should also be present.  If a batch effect is to be included in the
#' design, then this should be identified under the column <batch>.  If IP controls
#' are used then these should be specified in columns <input1> (and <input2> if paired
#' reads); this is not necessary for ATAC-seq.
#' @param reference character vector specifying the conditions in order.  For example,
#' c("A", "B", "C", "D") would mean "A" is the reference condition to which "B", "C"
#' and "D" are compared; in addition, "C" and "D" will be compared to "B", and "D"
#' will be compared to "C". If \code{NULL} then the comparisons will be arranged
#' alphabetically.  [DEFAULT = NULL].
#' @param species character string specifying the name of the species. Only
#' \code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].
#' @param output.dir character string specifying the directory to which results
#' will be saved.  If the directory does not exist, it will be created.
#' @param threads an integer value indicating the number of parallel threads to
#' be used by FastQC. [DEFAULT = maximum number of available threads - 1].
#'
#' @importFrom parallel detectCores
#'
#' @export
#'
#' @return A \code{data.frame}

sanity_check <- function(sample.info, reference = NULL,
                              species = c("human", "mouse"),
                              output.dir, threads = NULL) {
  inputdata <- list()

  # Results directory
  if (!dir.exists(output.dir)) create_dir(output.dir)
  inputdata$outdir <- output.dir

  # Threads
  if (is.null(threads)) threads <- parallel::detectCores() - 1
  inputdata$threads <- threads

  # Annotation
  species <- match.arg(species, c("human", "mouse"))
  inputdata$annotation <- species

  # Check sample info
  sd <- check_sample(sample.info, reference = reference)

  inputdata$sampleinfo <- sd$sample.info
  inputdata$input <- sd$input

  inputdata$paired <- sd$paired

  # Make design
  inputdata$design <- make_design(inputdata$sampleinfo)

  return(inputdata)
}
