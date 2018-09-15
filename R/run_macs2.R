#' Perform peakcalling with MACS2
#'
#' @param macs2 character string specifying path to macs2 [DEFAULT = "macs2"]
#' @param out.dir character string specifying the directory to which results
#' will be saved.  If the directory does not exist, it will be created. [DEFAULT = "peaks"]
#' @param treatment.files character vector specifying path to alignment files. [DEFAULT = NULL]
#' @param control.files character vector specifying path to input controls.  This is optional,
#' but if used then MUST be same length as <treatment.files>, so if the same input control is used
#' for multiple ltreatments, then it should be specified multiple times.  If NULL, then no.lambda
#' will be set to \code{TRUE}.  [DEFAULT = NULL]
#' @param input.format character vector specifying format of input file; can be "ELAND", "BED",
#' "ELANDMULTI", "ELANDEXPORT", "ELANDMULTIPET" (for pair-end tags), "SAM", "BAM", "BOWTIE",
#' "BAMPE" or "BEDPE". [DEFAULT = "AUTO"]
#' @param broad logical, if TRUE then MACS2 will output broadpeak results
#' @param species character string specifying the name of the species. Only
#' \code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].
#' @param sample.names character vector specifying sample names to use as prefix for files. If NULL,
#' this will be derived from the file names.  [DEFAULT = NULL]
#' @param no.model logical, if TRUE then MACS2 shifting model will be bypassed.
#' @param extsize numeric vector.  If no.model = TRUE, this parameter is to extend reads in 5'->3'
#' direction to fix-sized fragments.  Ignored if no.model = FALSE or paired end reads. If used, MUST
#' be same length as <treatment.files> and in same order. [DEFAULT = 0]
#' @param shift numeric vector. If no.model = TRUE, this value is used to move cutting ends (5') then
#' apply extsize from 5' to 3' direction to extend them to fragments. When this value is negative,
#' ends will be moved toward 3'->5' direction, otherwise 5'->3' direction. [DEFAULT = 0]
#' @param no.lambda logical, if TRUE then the background lambda will be used as local lambda.
#' @param call.summits logical, if TRUE then MACS2 will reanalyze the shape of signal profile to deconvolve
#' subpeaks within each peak. If used, the output subpeaks of a big peak region will have the same peak boundaries,
#' but different scores and peak summit positions. [DEFAULT = TRUE]
#' @param qvalue numeric value specifying the minimum FDR cutoff to call significant regions. [DEFAULT = 0.05]
#' @param verbosity integer value specifying verbosity.  If 0 then only critical messages displayed.  If greater
#' or equal to 3 then all messages shown.
#' @param threads an integer value indicating the number of workers to be used. If NULL then one less than the
#' maximum number of cores will be used. [DEFAULT = NULL].
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%

run_macs2 <- function(macs2 = "macs2", out.dir = "peaks",
                      treatment.files = NULL, control.files = NULL, input.format = "AUTO", broad = FALSE,
                      species = c("human", "mouse"), sample.names = NULL, no.model = FALSE, extsize = 0, shift = 0,
                      no.lambda = FALSE, call.summits = TRUE, qvalue = 0.05, verbosity = 3, threads = NULL) {

  i <- NULL

  species <- match.arg(species)
  if (species == "human") genome.size <- "hs"
  if (species == "mouse") genome.size <- "mm"

  create_dir(out.dir)

  # create output name string
  if (is.null(sample.names)) {
    sample.names <- reduce_path(treatment.files)
  }
  sample.names2 <- file.path(out.dir, paste0(sample.names, "_fdr", 100 * qvalue))

  cmd.string <- sprintf("%s callpeak -t %s -f %s -g %s -n %s -B -q %s --verbose %s",
                        macs2, treatment.files, input.format, genome.size, sample.names2, qvalue, verbosity)

  if (is.null(control.files)) no.lambda <- TRUE
  if (no.lambda) cmd.string <- paste(cmd.string, "--nolambda")

  if (!is.null(control.files)) {
    if (length(control.files) != length(treatment.files)) {
      give_error("/nThe number of control and treatment files do not match./n")
    } else {
      cmd.string <- paste(cmd.string, "-c", control.files)
    }
  }

  if (no.model) {
    cmd.string <- paste(cmd.string, "--nomodel --extsize", extsize, "--shift", shift)
  }

  if (call.summits) cmd.string <- paste(cmd.string, "--call-summits")
  if (broad) cmd.string <- paste(cmd.string, "--broad")

  give_note("Running MACS2...\n\n")

  if (is.null(threads)) threads <- parallel::detectCores() - 1

  if (threads >= 2) {
    cl <- parallel::makeCluster(threads)
    doParallel::registerDoParallel(cl)
    foreach::foreach(i = seq_along(cmd.string)) %dopar% {
      run_cmd <- function(cmd, intern = FALSE) {
        if (.Platform$OS.type != "windows") {
          system(command = cmd, intern = intern)
        } else {
          shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
        }
      }
      cat(cmd.string)
      run_cmd(cmd.string[[i]])
    }
    parallel::stopCluster(cl)
  } else {
    foreach::foreach(i = seq_along(cmd.string)) %do% {
      run_cmd <- function(cmd, intern = FALSE) {
        if (.Platform$OS.type != "windows") {
          system(command = cmd, intern = intern)
        } else {
          shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
        }
      }
      cat(cmd.string)
      run_cmd(cmd.string[[i]])
    }
  }

  peak.files <- list.files(out.dir, "broadPeak|narrowPeak", full.names = TRUE)
  peaks <- lapply(peak.files, peak2Granges)
  names(peaks) <- sample.names
  return(peaks)
}
