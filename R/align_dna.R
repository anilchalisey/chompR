#' Alignment of DNA-seq reads
#'
#' @inheritParams run_chip
#' @param reads1 Character vector of mate1 reads.  If specified, then reads.dir
#'   must be NULL.
#' @param reads2 Character vector of mate2 reads.  If specified, then reads.dir
#'   must be NULL.  Must be the same length as mate1.  If single-end sequencing,
#'   then should be left as NULL.
#' @param fastq Logical indicating if reads are FASTQ files.
#' @param fasta Logical indicating if reads are FASTA files.
#' @param softClipPenalty Sets the maximum (MX) and minimum (MN) penalties for
#' soft-clipping per base, both integers.  Must be given in the format "MX,MN".
#' @param noSoftClip Logical indicating whether to disallow soft-clipping.
#' @param tmo Logical indicating whether to report only those reads aligning to
#' known transcripts.
#' @param maxAlign Integer indicating the maximum number of distinct primary
#' alignments to search for each read.
#' @param secondary Logical indicating whether to report secondary alignments.
#' @param nomixed By default, when hisat2 cannot find a concordant or discordant
#' alignment for a pair, it then tries to find alignments for the individual
#' mates. If TRUE, this option disables that behavior.
#' @param nodiscordant By default, hisat2 looks for discordant alignments if it
#' cannot find any concordant alignments.  If true, this option disables that
#' behavior.
#' @param rgid Character string, to which the read group ID is set.
#' @param quiet If TRUE, print nothing except alignments and serious errors.
#' @param non_deterministic When set to TRUE, HISAT2 re-initializes its
#' pseudo-random generator for each read using the current time.
#' @param maxInsert The maximum fragment length for valid paired-end alignments.
#' This option is valid only with noSplice = TRUE.
#' @param memory String specifying maximum memory per thread; suffix K/M/G
#'   recognized.
#' @param remove.mitochondrial Character string.  If set, this will remove reads
#'   mapping to the mitochondrial genome.  The string should match the reference
#'   name for the mitochindrial genome in the alignment file.  Examples include
#'   "ChrM", "M" and "MT".
#' @param remove.duplicates If TRUE, duplicate reads will be removed.
#' @param hash_table Size of hash table for finding read pairs (default is
#' 262144 reads); will be rounded down to the nearest power of two.  For best
#' performance should be > (average coverage) * (insert size).
#' @param overflow_size Size of the overflow list where reads, thrown out of
#' the hash table, get a second chance to meet their pairs (default is
#' 200000 reads); increasing the size reduces the number of temporary files
#' created.
#' @param io_buffer Controls sizes of the two buffers (in MB) used for reading
#' and writing BAM during the second pass (default is 128).
#'
#' @return Raw and filtered BAM files
#'
#' @export

align_dna <- function(
  ## Important - general options
  threads = 1,
  output.dir = ".",
  hisat2 = "hisat2",
  samtools = "samtools",
  sambamba = "sambamba",
  species = c("human", "mouse"),
  ## Important - align functions
  idx = NULL,
  reads1 = NULL,
  reads2 = NULL,
  ## Can usually be left as default - align functions
  fastq = TRUE,
  fasta = FALSE,
  softClipPenalty = NULL,
  noSoftClip = FALSE,
  tmo = FALSE,
  secondary = FALSE,
  maxAlign = NULL,
  nomixed = FALSE,
  nodiscordant = FALSE,
  rgid = NULL,
  quiet = FALSE,
  non_deterministic = FALSE,
  maxInsert = NULL,
  # Can usually be left as default - samtools functions
  memory = "1G",
  remove.mitochondrial = "MT",
  remove.duplicates = TRUE,
  # Can usually be left as default - sambamba function
  hash_table = 262144,
  overflow_size = 200000,
  io_buffer = 128) {

  if (check_cmd(hisat2) == "Not Found") {
    give_error("the hisat2 path is invalid or it is not installed correctly.")
  }

  if (check_cmd(samtools) == "Not Found") {
    give_error("the samtools path is invalid or it is not installed correctly.")
  }

  if (check_cmd(sambamba) == "Not Found") {
    give_error("the sambamba path is invalid or it is not installed correctly.")
  }

  # Get the file names
  mate1 <- reads1
  mate2 <- reads2

  if (!is.null(mate2)) {
    if (length(mate1) != length(mate2)) {
      give_error("The number of reads1 files and reads2 files do not match.")
    }
  }

  species <- match.arg(species)

  if (is.null(idx)) {
    idx <- build_index(species = species)
  }


  # Align
  sam.files <-
    run_hisat2(hisat2 = hisat2, idx = idx, mate1 = mate1, mate2 = mate2,
               fastq = fastq, fasta = fasta, softClipPenalty = softClipPenalty,
               noSoftClip = noSoftClip, noSplice = TRUE,
               knownSplice = NULL, strand = NULL, tmo = tmo,
               maxAlign = maxAlign, secondary = secondary, minInsert = NULL,
               maxInsert = maxInsert, nomixed = nomixed,
               nodiscordant = nodiscordant, threads = threads, rgid = rgid,
               quiet = quiet, non_deterministic = non_deterministic)

  # Convert to BAM
  bam.files <-
    run_samsort(samtools = samtools, file = sam.files, outformat = "BAM",
                threads = threads, memory = memory, sortbyname = FALSE,
                suffix = "", keep = FALSE)

  # Mark duplicates and index
  run_sambambadup(sambamba = sambamba, bamfile = bam.files,
                  outfile = NULL, remove = FALSE,
                  threads = threads, hash_table = hash_table,
                  overflow_size = overflow_size, io_buffer = io_buffer)

  unlink(bam.files)
  torename <- paste0(reduce_path(bam.files), "_markdup.bam")
  file.rename(from = torename, to = bam.files)
  file.rename(from = paste0(torename, ".bai"), to = paste0(bam.files, ".bai"))

  # Get some statistics for the BAM files
  give_note("\n\nGenerating alignment stats from BAM files\n\n")
  diag.stats <-
    run_samflagstat(samtools = samtools, threads = threads, bamfile = bam.files)
  cnames <- c("", basename(colnames(diag.stats)))
  utils::write.table(diag.stats, file = "alignment_stats.txt",
                     col.names = NA,
                     row.names = TRUE, quote = FALSE, sep = "\t")

  # Filter out mitochondrial and improperly aligned reads

  if (!is.null(mate2)) {
    filtered.bam.files <-
      run_samview(samtools = samtools, file = bam.files, regions = NULL,
                  chrom.sizes = NULL, include.flag = NULL, exclude.flag = NULL,
                  minQual = NULL, outformat = "BAM", outname = NULL,
                  include.header = FALSE, count = FALSE, threads = threads,
                  subsample = NULL, keep.paired = TRUE, keep.proper.pair = TRUE,
                  remove.unmapped = TRUE, remove.not.primary = TRUE,
                  remove.duplicates = remove.duplicates,
                  remove.supplementary.alignment = TRUE,
                  remove.mitochondrial = shQuote(paste0("KI\\|GL\\|\\*\\|", remove.mitochondrial)))
  } else {
    filtered.bam.files <-
      run_samview(samtools = samtools, file = bam.files, regions = NULL,
                  chrom.sizes = NULL, include.flag = NULL, exclude.flag = NULL,
                  minQual = NULL, outformat = "BAM", outname = NULL,
                  include.header = FALSE, count = FALSE, threads = threads,
                  subsample = NULL, keep.paired = FALSE,
                  keep.proper.pair = FALSE, remove.unmapped = TRUE,
                  remove.not.primary = TRUE,
                  remove.duplicates = remove.duplicates,
                  remove.supplementary.alignment = TRUE,
                  remove.mitochondrial = shQuote(paste0("KI\\|GL\\|\\*\\|", remove.mitochondrial)))
  }

  # create index
  run_samindex(samtools = samtools, bamfile = filtered.bam.files,
               threads = threads)

  bamdir <- file.path(output.dir, "bam")
  dir.create(bamdir, recursive = TRUE, showWarnings = FALSE)
  lapply(bam.files, function(x) file.rename(x, file.path(bamdir, x)))
  bam.bai <- paste0(bam.files, ".bai")
  lapply(bam.bai, function(x) file.rename(x, file.path(bamdir, x)))
  lapply(list.files(pattern = "*.log"),
         function(x) file.rename(x, file.path(bamdir, x)))

  filtereddir <- file.path(output.dir, "filteredbam")
  dir.create(filtereddir, recursive = TRUE, showWarnings = FALSE)
  filteredbam <- list.files(pattern = "filtered.bam")
  lapply(filteredbam, function(x) {
    file.rename(x, file.path(filtereddir, gsub("_filtered", "", x)))
  })

  file.rename(from = "alignment_stats.txt",
              to = file.path(output.dir, "alignment_stats.txt"))

  return(invisible())
}

