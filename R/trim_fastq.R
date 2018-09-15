#' An R-based wrapper for fastp
#'
#' @description Run the fastp tool
#'
#' @section fastp path:
#' If the executable is in \code{$PATH}, then the default value for paths
#' (\code{"fastp"}) will work. If it is not in \code{$PATH}, then the absolute
#' path should be given. If using Windows 10, it is assumed that fastp has
#' been installed in WSL, and the same rules apply.
#'
#'
#' @details This script runs the fastp tool and requires installation
#' of fastp.  Pre-compiled binaries and installation instructions may
#' be found at \url{https://github.com/OpenGene/fastp}
#'
#' @references
#' Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu (2018):
#' fastp: an ultra-fast all-in-one FASTQ preprocessor.
#' BioRxiv 274100; \url{https://doi.org/10.1101/274100}
#'
#' @param fastp a character string specifying the path to the fastp executable.
#' [DEFAULT = "fastp"].
#' @param fastq1 a character vector indicating the read files to be trimmed.
#' @param fastq2 (optional) a character vector indicating read files to be
#' trimmmed.  If specified, it is assumed the reads are paired, and this vector
#' MUST be in the same order as those listed in \code{fastq1}.  If \code{NULL}
#' then it is assumed the reads are single-end. [DEFAULT = NULL]
#' @param dest.dir a character string specifying the output directory.  If NULL
#' a directory named "TRIMMED_FASTQC" is created in the current working directory.
#' [DEFAULT = NULL].
#' @param disable_adapter_trimming logical, if TRUE adapter trimming is disabled. [DEFAULT = FALSE]
#' @param adapter_sequence character string, specifying the adapter for read1. For SE data, if not specified,
#' the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. [DEFAULT = NULL]
#' @param adapter_sequence_r2 character string, the adapter for read2 (PE data only). This is used if R1/R2 are
#' found not overlapped. If not specified, it will be the same as <adapter_sequence>. [DEFAULT = NULL]
#' @param trim_front1 integer specifying number of bases to trim at 5' end for read1 [DEFAULT = 0]
#' @param trim_front2 integer specifying number of bases to trim at 5' end for read2 [DEFAULT = 0]
#' @param trim_tail1 integer specifying number of bases to trim at 3' end for read1 [DEFAULT = 0]
#' @param trim_tail2 integer specifying number of bases to trim at 3' end for read2 [DEFAULT = 0]
#' @param trim_poly_g logical, if TRUE, force polyG tail trimming [DEFAULT = FALSE].
#' @param poly_g_min_len integer specifying the minimum length to detect polyG in the read tail. [DEFAULT = 10]
#' @param trim_poly_x logical, f TRUE, enable polyX trimming in 3' ends. [DEFAULT = FALSE]
#' @param poly_x_min_len integer specifying the minimum length to detect polyX in the read tail. [DEFAULT = 10]
#' @param cut_by_quality5 logical, if TRUE enable per read cutting by quality at 5' end (WARNING: this will
#' interfere deduplication for both PE/SE data) [DEFAULT = FALSE]
#' @param cut_by_quality3 logical, if TRUE enable per read cutting by quality at 3' end (WARNING: this will
#' interfere deduplication for both SE data) [DEFAULT = FALSE]
#' @param cut_window_size integer specifying the base pair size of the sliding window for sliding window trimming [DEFAULT = 4]
#' @param cut_mean_quality integer specifying the mean phred quality threshold within a sliding window for removing
#' bases [DEFAULT = 20]
#' @param disable_quality_filtering logical, if TRUE then quality filtering is enabled. [DEFAULT = TRUE]
#' @param qualified_quality_phred integer specifying the base quality threshold. [DEFAULT = 15]
#' @param unqualified_percent_limit numeric specifying the percentage of bases allowed to be below the threshold
#' before a read/pair is discarded. [DEFAULT = 40]
#' @param n_base_limit integer specifying the number of allowable uncallable reads (N) before a read/pair is discarded.
#' [DEFAULT = 5]
#' @param disable_length_filtering logical, if TRUE then length filtering is enabled. [DEFAULT = TRUE]
#' @param length_required integer specifying the length below which reads will be discarded. [DEFAULT = 15]
#' @param length_limit integer specifying the length above which reads will be discarded; if 0 then no limit applied. [DEFAULT = 0]
#' @param low_complexity_filter logical, if TRUE then enable low complexity filter. The complexity is defined as the percentage of
#' base that is different from its next base (base[i] != base[i+1]). [DEFAULT = FALSE]
#' @param complexity_threshold numeric specifying the threshold for the low complexity filter (0~100). [DEFAULT = 30]
#' @param filter_by_index1  character string specifying a file containing a list of barcodes of index1 to be filtered
#' out, one barcode per line. [DEFAULT = NULL]
#' @param filter_by_index2 character string specifying a file containing a list of barcodes of index2 to be filtered
#' out, one barcode per line. [DEFAULT = NULL]
#' @param filter_by_index_threshold the allowed difference of index barcode for index filtering; 0 means completely
#' identical. [DEFAULT = 0]
#' @param correction logical, if TRUE theb enable base correction in overlapped regions (only for PE data). [DEFAULT = FALSE]
#' @param overlap_len_require integer specifying the minimum length of the overlapped region for overlap analysis
#' based adapter trimming and correction. [DEFAULT = 30]
#' @param overlap_diff_limit integer specifying the maximum difference of the overlapped region for overlap analysis
#' based adapter trimming and correction. [DEFAULT = 5]
#' @param overrepresentation_analysis logical, if TRUE then enable overrepresented sequence analysis. [DEFAULT = FALSE]
#' @param overrepresentation_sampling integer specifying how reads will be computed for overrepresentation analysis, e.g.
#' if set to 20, then 1-in029 reads will be sampled.  May range from 1 to 10000; smaller is slower, [DEFAULT = 20]
#' @param threads an integer value indicating the number of workers to be used. If NULL then one less than the
#' maximum number of cores will be used. [DEFAULT = NULL].
#'
#' @importFrom parallel detectCores
#'
#' @export

trim_fastq <- function(fastp = "fastp", fastq1, fastq2 = NULL, dest.dir = NULL,
                        disable_adapter_trimming = FALSE, adapter_sequence = NULL,
                        adapter_sequence_r2 = NULL, trim_front1 = 0, trim_front2 = 0,
                        trim_tail1 = 0, trim_tail2 = 0, trim_poly_g = FALSE,
                        poly_g_min_len = 10, trim_poly_x = FALSE, poly_x_min_len = 10,
                        cut_by_quality5 = FALSE, cut_by_quality3 = FALSE,
                        cut_window_size = 4, cut_mean_quality = 20,
                        disable_quality_filtering = FALSE, qualified_quality_phred = 15,
                        unqualified_percent_limit = 40, n_base_limit = 5,
                        disable_length_filtering = FALSE, length_required = 15,
                        length_limit = 0, low_complexity_filter = FALSE,
                        complexity_threshold = 30, filter_by_index1 = NULL,
                        filter_by_index2 = NULL, filter_by_index_threshold = 0,
                        correction = FALSE, overlap_len_require = 30, overlap_diff_limit = 5,
                        overrepresentation_analysis = FALSE, overrepresentation_sampling = 20,
                        threads = NULL) {

  if (is.null(fastq2)) {
    paired <- FALSE
  } else {
    if(length(fastq1) != length(fastq2)) {give_error("The number of forward and reverse reads do not match")}
    paired <- TRUE
  }

  if (is.null(dest.dir)) dest.dir <- "TRIMMED_FASTQC"
  dir.create(dest.dir, showWarnings = FALSE)
  cmd <- fastp
  output1 <- file.path(dest.dir, basename(fastq1))
  output2 <- file.path(dest.dir, basename(fastq2))

  if (paired) {
    give_note("\nPaired end reads detected...\n\n")
    if (disable_adapter_trimming) cmd <- paste(cmd, "-A")
    if (!is.null(adapter_sequence)) cmd <- paste(cmd, "--adapter_sequence", adapter_sequence)
    if (!is.null(adapter_sequence_r2)) cmd <- paste(cmd, "--adapter_sequence_r2", adapter_sequence_r2)
    cmd <- sprintf("%s -f %s -t %s -F %s -T %s", cmd, trim_front1, trim_tail1, trim_front2, trim_tail2)
    if (!trim_poly_g) cmd <- paste(cmd, "-g")
    if (trim_poly_g) cmd <- paste(cmd, "poly_g_min_len", poly_g_min_len)
    if (trim_poly_x) cmd <- paste(cmd, "-x --poly_x_min_len", poly_x_min_len)
    if (cut_by_quality5) cmd <- paste(cmd, "-5")
    if (cut_by_quality3) cmd <- paste(cmd, "-3")
    if (cut_by_quality5 | cut_by_quality3) cmd <- paste(cmd, "-W", cut_window_size, "-M", cut_mean_quality)
    if (disable_quality_filtering) {
      cmd <- paste(cmd, "-Q")
    } else {
      cmd <- paste(cmd, "--qualified_quality_phred", qualified_quality_phred,
                   "--unqualified_percent_limit", unqualified_percent_limit,
                   "--n_base_limit", n_base_limit)
    }
    if (disable_length_filtering) {
      cmd <- paste(cmd, "-L")
    } else {
      cmd <- paste(cmd, "--length_required", length_required,
                   "--length_limit", length_limit)
    }
    if (low_complexity_filter) cmd <- paste(cmd, "-y -Y", complexity_threshold)
    if (!is.null(filter_by_index1)) cmd <- paste(cmd, "--filter_by_index1")
    if (!is.null(filter_by_index2)) cmd <- paste(cmd, "--filter_by_index2")
    if (!is.null(filter_by_index1) | !is.null(filter_by_index1)) {
      cmd <- paste(cmd, "--filter_by_index_threshold", filter_by_index_threshold)
    }
    if (correction) cmd <- sprintf("%s -c --overlap_len_require %s --overlap_diff_limit %s",
                                   cmd, overlap_len_require, overlap_diff_limit)
    if (overrepresentation_analysis) cmd <- paste(cmd, "-p -P", overrepresentation_sampling)
    if (is.null(threads)) threads <- parallel::detectCores() - 1
    cmd <- paste(cmd, "-w", threads)

  } else {

    if (disable_adapter_trimming) cmd <- paste(cmd, "-A")
    if (!is.null(adapter_sequence)) cmd <- paste(cmd, "--adapter_sequence", adapter_sequence)
    cmd <- sprintf("%s -f %s -t %s", cmd, trim_front1, trim_tail1)
    if (!trim_poly_g) cmd <- paste(cmd, "-g")
    if (trim_poly_g) cmd <- paste(cmd, "poly_g_min_len", poly_g_min_len)
    if (trim_poly_x) cmd <- paste(cmd, "-x --poly_x_min_len", poly_x_min_len)
    if (cut_by_quality5) cmd <- paste(cmd, "-5")
    if (cut_by_quality3) cmd <- paste(cmd, "-3")
    if (cut_by_quality5 | cut_by_quality3) cmd <- paste(cmd, "-W", cut_window_size, "-M", cut_mean_quality)
    if (disable_quality_filtering) {
      cmd <- paste(cmd, "-Q")
    } else {
      cmd <- paste(cmd, "--qualified_quality_phred", qualified_quality_phred,
                   "--unqualified_percent_limit", unqualified_percent_limit,
                   "--n_base_limit", n_base_limit)
    }
    if (disable_length_filtering) {
      cmd <- paste(cmd, "-L")
    } else {
      cmd <- paste(cmd, "--length_required", length_required,
                   "--length_limit", length_limit)
    }
    if (low_complexity_filter) cmd <- paste(cmd, "-y -Y", complexity_threshold)
    if (!is.null(filter_by_index1)) cmd <- paste(cmd, "--filter_by_index1")
    if (!is.null(filter_by_index2)) cmd <- paste(cmd, "--filter_by_index2")
    if (!is.null(filter_by_index1) | !is.null(filter_by_index1)) {
      cmd <- paste(cmd, "--filter_by_index_threshold", filter_by_index_threshold)
    }
    if (overrepresentation_analysis) cmd <- paste(cmd, "-p -P", overrepresentation_sampling)
    if (is.null(threads)) threads <- parallel::detectCores() - 1
    cmd <- paste(cmd, "-w", threads)

  }

  give_note("\nRemoving adapters and performing quality trimming...\n\n")

  result <- vector("list", length(fastq1))

  if (paired) {
    for (i in seq_along(fastq1)) {
      cmd <- sprintf("%s -i %s -I %s -o %s -O %s", cmd, fastq1[[i]], fastq2[[i]],
                     output1[[i]], output2[[i]])
      result[[i]] <- run_cmd(cmd, intern = TRUE)
      run_cmd(sprintf("mv fastp.html %s",
                      file.path(dest.dir, paste0(reduce_path(fastq1[[i]]), "-", reduce_path(fastq2[[i]]), ".fastp.html"))))
      run_cmd(sprintf("mv fastp.json %s",
                      file.path(dest.dir, paste0(reduce_path(fastq1[[i]]), "-", reduce_path(fastq2[[i]]), ".fastp.json"))))
    }
  } else {
      for (i in seq_along(fastq1)) {
        cmd <- sprintf("%s -i %s -o %s", cmd, fastq1[[i]], output1[[i]])
        result[[i]] <- run_cmd(cmd, intern = TRUE)
        run_cmd(sprintf("mv fastp.html %s",
                        file.path(dest.dir, paste0(reduce_path(fastq1[[i]]), "-", reduce_path(fastq2[[i]]), ".fastp.html"))))
        run_cmd(sprintf("mv fastp.json %s",
                        file.path(dest.dir, paste0(reduce_path(fastq1[[i]]), "-", reduce_path(fastq2[[i]]), ".fastp.json"))))
      }
  }

  for (i in seq_along(fastq1)) {
    if (paired) {
      names(result)[[i]] <- paste0(reduce_path(fastq1[[i]]), "-", reduce_path(fastq2[[i]]))
      } else {
        names(result)[[i]] <- paste0(reduce_path(fastq1[[i]]))
      }
  }
  result <- lapply(names(result), function(x) c(paste("\n\n*****\n", x), result[[x]]))
  lapply(result, write, file.path(dest.dir, "fastp.log"), append = TRUE)
  return(invisible())
}




