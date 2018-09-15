#' Check that targets file contains the necessary info in the right format
#'
#' @inheritParams sanity_check
#'
#' @importFrom utils read.table
#' @importFrom stats relevel

check_sample <- function(sample.info, reference = NULL) {
  if (!file.exists(sample.info)) {
    give_error("SAMPLE.INFO ERROR:\nsample information file does not exist or path is incorrect.")
  }
  sample.info <- utils::read.table(sample.info, header = TRUE, stringsAsFactors = FALSE)
  names(sample.info) <- tolower(names(sample.info))

  for (required in c("condition", "sample", "file1")) {
    if (!required %in% names(sample.info)) {
      give_error(paste0("SAMPLE.INFO ERROR:\nCould not find required field '",
                         required, "' in sample information file."))
    }
  }

  metadata <- sample.info

  input <- ifelse("input1" %in% names(metadata), TRUE, FALSE)
  if (input) give_note("Input/control files detected")

  if (min(table(metadata$condition)) < 2) {
    give_error("SAMPLE.INFO ERROR:\nLess than 2 replicates in smallest group from <condition>.")
  }
  if (length(unique(metadata$condition)) < 2) {
    give_error("SAMPLE.INFO ERROR:\nfield <condition> needs at least two different values/groups.")
  }

  for (unique_values_required_field in c("sample", "file1")) {
    if (anyDuplicated(sample.info[, unique_values_required_field])) {
      give_error(paste("SAMPLE.INFO ERROR:\nValues in field",
                       unique_values_required_field,
                       "in sample info file are not unique."))}
  }

  for (path in sample.info$file1) {
    if (!file.exists(path)) {
      give_error(paste("\nIncorrect path to", path, "\nFile not found.\n"))
    }
  }

  if (!all(grepl("fastq|fq", sample.info$file1))) {
    give_error("ERROR:\n Paths in sample info file should all be <.fastq/.fq>\n")
  }

  if ("file2" %in% names(sample.info)) {
    if (!all(grepl("fastq|fq", sample.info$file2))) {
        give_error("ERROR:\n Paths in sample info file should all be <.fastq/.fq>\n")
    }

    for (path in sample.info$file2) {
        if (!file.exists(path)) {
          give_error(paste0("\nIncorrect path to ", path, ". File not found.\n"))
        }
      }

    give_note("\nEach sample has 2 files - assuming PE reads\n")
    paired <- TRUE
  } else {
    give_note("\nOnly one file per sample - assuming SE reads\n")
    paired <- FALSE
  }

  if ("input1" %in% names(sample.info)) {
    if (!all(grepl("fastq|fq", sample.info$input1))) {
      give_error("ERROR:\n Paths in sample info file should all be <.fastq/.fq>\n")
    }

    for (path in sample.info$input1) {
      if (!file.exists(path)) {
        give_error(paste0("\nIncorrect path to ", path, ". File not found.\n"))
      }
    }
  }

  if ("input2" %in% names(sample.info)) {
    if (!all(grepl("fastq|fq", sample.info$input2))) {
      give_error("ERROR:\n Paths in sample info file should all be <.fastq/.fq>\n")
    }

    for (path in sample.info$input2) {
      if (!file.exists(path)) {
        give_error(paste0("\nIncorrect path to ", path, ". File not found.\n"))
      }
    }
  }

  if (is.null(reference)) {
    metadata$condition <- factor(metadata$condition)
  } else {
      if (length(reference) == 1) {
      metadata$condition <- stats::relevel(factor(metadata$condition), ref = reference)
    } else {
      if (length(reference) > 1) {
        if (length(setdiff(metadata$condition, reference)) > 0) {
          give_error("ERROR:\n The factors specified in <reference> are those in the condition column do not match\n")
        }
        if (length(setdiff(reference, metadata$condition)) > 0) {
          give_error("ERROR:\n The factors specified in <reference> are those in the condition column do not match\n")
        }
        metadata$condition <- factor(metadata$condition, levels = reference)
      }
    }
  }

  if ("batch" %in% names(metadata)) metadata$batch <- as.factor(metadata$batch)

  return(list(sample.info = metadata,
              paired = paired, input = input))

}
