#' Run ChIP-seq analysis
#'
#' Complete pipeline from quality assessment of raw reads through to peak calling
#' and differential analysis.
#'
#' @param sample.info character string giving the path to a tab-delimited text
#' file with at least the columns <condition> (treatment condition),
#' <sample> (sample name), and <file1> (absolute or relative path to the
#' fastq files).  If fastq files and PE reads, then a column
#' <file2> should also be present.  If a batch effect is to be included in the
#' design, then this should be identified under the column <batch>.  If IP controls
#' are used then these should be specified in columns <input1> (and <input2> if paired
#' reads).
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
#' @param bigwig logical, if \code{TRUE} then bigwif files will be created.
#' @param fastqc a character string specifying the path to the fastqc executable.
#' [DEFAULT = "fastqc"].
#' @param multiqc a character string specifying the path to the multiqc executable.
#' [DEFAULT = "multiqc"].
#' @param fastp a character string specifying the path to the fastp executable.
#' [DEFAULT = "fastp"].
#' @param samtools a character string specifying the path to the samtools executable.
#' [DEFAULT = "samtools"].
#' @param hisat2 a character string specifying the path to the hisat2 executable.
#' [DEFAULT = "hisat2"].
#' @param sambamba a character string specifying the path to the sambamba executable.
#' [DEFAULT = "sambamba"].
#' @param macs2 a character string specifying the path to the MACS2 executable.
#' [DEFAULT = "macs2"].
#' @param featurecounts a character string specifying the path to the featurecounts executable.
#' [DEFAULT = "featureCounts"].
#' @param idx character vector specifying the basename of the index for the reference genome.
#' The basename is the name of any of the index files up to but not including the final
#' .1.ht2, etc. If \code{NULL} then the index for the relevant species (human or mouse) will be
#' created using the \code{build_index()} function.
#' @param broad logical, if \code{TRUE} then broad peaks will be called.
#'
#' @export

run_chip <- function(# Important options
  sample.info,
  reference = NULL,
  species = c("human", "mouse"),
  output.dir,
  threads = NULL,
  bigwig = FALSE,
  # QC options
  fastqc = "fastqc",
  multiqc = "multiqc",
  fastp = "fastp",
  # Alignment options
  samtools = "samtools",
  hisat2 = "hisat2",
  sambamba = "sambamba",
  idx = NULL,
  # Other options
  macs2 = "macs2",
  featurecounts = "featureCounts",
  broad = FALSE) {

  file1 <- . <- input1 <- bam.x <- filteredbam.x <- NULL
  bam.y <- filteredbam.y <- filteredbam <- NULL
  Var2 <- value <- Var1 <- peak <- NULL

  ## STEP 1: Import the metadata information ----------------------------------
  give_note("\nImporting metadata\n")
  metadata <-
    sanity_check(
      sample.info = sample.info,
      reference = reference,
      species = species,
      output.dir = output.dir,
      threads = threads
    )


  ## STEP 2: QC of reads ------------------------------------------------------
  # Run FastQC
  # Code is written so that if the files already exist, then the function is skipped
  fqc.out <-
    c(
      metadata$sampleinfo$file1,
      metadata$sampleinfo$file2,
      unique(metadata$sampleinfo$input1),
      unique(metadata$sampleinfo$input2)
    )
  fqc.out.path <-
    paste0(file.path(metadata$outdir, "fastqc", reduce_path(fqc.out)),
           "_fastqc.zip")
  if (!all(file.exists(fqc.out.path))) {
    give_note("Running QC of FASTQ files\n")
    run_fastqc(
      fastqc.files = fqc.out,
      dest.dir = file.path(metadata$outdir, "fastqc"),
      threads = metadata$threads,
      fastqc = fastqc
    )
  }

  # Run multiQC
  multiqc.out.path <-
    file.path(metadata$outdir, "multiqc", "multiqc_fastqc.html")
  if (!file.exists(multiqc.out.path)) {
    run_multiqc(
      fastqc.dir = file.path(metadata$outdir, "fastqc"),
      dest.dir = file.path(metadata$outdir, "multiqc"),
      multiqc = multiqc
    )
  }

  # Trim and filter reads
  trimfq.out.path <-
    paste0(file.path(metadata$outdir, "trimmed", reduce_path(fqc.out)),
           ".fastq")
  if (!all(file.exists(trimfq.out.path))) {
    give_note("Trimming and filtering reads\n")
    trim_fastq(
      fastq1 = c(
        metadata$sampleinfo$file1,
        unique(metadata$sampleinfo$input1)
      ),
      fastq2 = c(
        metadata$sampleinfo$file2,
        unique(metadata$sampleinfo$input2)
      ),
      dest.dir = file.path(metadata$outdir, "trimmed"),
      fastp = fastp
    )
  }

  ## STEP 3: Align trimmed/filtered reads.  Remove duplicated reads and -------
  ## those mapping to mitochondrial and non-standard chromosomes --------------
  # Use trimmed reads for alignment
  bam.out <-
    c(metadata$sampleinfo$file1,
      unique(metadata$sampleinfo$input1)) %>%
    reduce_path() %>% sub("_1$", "", .) %>% paste0(".bam")

  bam.out.path <-
    c(
      file.path(metadata$outdir, "bam", bam.out),
      file.path(metadata$outdir, "filteredbam", bam.out)
    )

  if (!all(file.exists(bam.out.path))) {
    give_note("Aligning reads\n")
    fq1 <- file.path(metadata$outdir, "trimmed",
                     basename(c(
                       metadata$sampleinfo$file1,
                       unique(metadata$sampleinfo$input1)
                     )))
    if (metadata$paired) {
      fq2 <- file.path(metadata$outdir, "trimmed",
                       basename(c(
                         metadata$sampleinfo$file2,
                         unique(metadata$sampleinfo$input2)
                       )))
      align_dna(
        threads = metadata$threads,
        output.dir = metadata$outdir,
        hisat2 = hisat2,
        samtools = samtools,
        sambamba = sambamba,
        species = metadata$annotation,
        idx = idx,
        reads1 = fq1,
        reads2 = fq2
      )
    } else {
      align_dna(
        threads = metadata$threads,
        output.dir = metadata$outdir,
        hisat2 = hisat2,
        samtools = samtools,
        sambamba = sambamba,
        species = metadata$annotation,
        idx = idx,
        reads1 = fq1,
        reads2 = fq2
      )
    }
  }

  index.out.path <- paste0(bam.out.path, ".bai")
  index.out.actual <-
    c(
      list.files(
        file.path(metadata$outdir, "bam"),
        pattern = "bai",
        full.names = TRUE
      ),
      list.files(
        file.path(metadata$outdir, "filteredbam"),
        pattern = "bai",
        full.names = TRUE
      )
    )

  if (length(setdiff(index.out.path, index.out.actual)) >= 1) {
    run_samindex(
      samtools = samtools,
      bamfile = bam.out.path,
      threads = metadata$threads
    )
  }

  ## STEP 4: QC of aligned reads ----------------------------------------------
  # Calculate some alignment statistics

  if (!file.exists(file.path(metadata$outdir, "bam", "diagstats.txt"))) {
    give_note("Performing QC of aligned reads\n")
    qc <-
      qc_alignment2(
        bam.files = list.files(file.path(metadata$outdir, "bam"), "*.bam$",
                               full.names = TRUE),
        threads = metadata$threads
      )
    run_cmd(sprintf("mv diagstats.txt %s", file.path(metadata$outdir, "bam")))
  } else {
    qc <-
      read.table(
        file.path(metadata$outdir, "bam", "diagstats.txt"),
        sep = "\t",
        header = TRUE
      )
  }

  # Perform cross-correlation on the (deduplicated) ChIP samples

  if (!file.exists(file.path(metadata$outdir, "cross-correlation", "length_results.txt"))) {
    give_note("Performing cross-correlation analysis\n")
    cc <-
      cc_chip(
        bam.files = list.files(file.path(metadata$outdir, "bam"),
                               "*.bam$", full.names = TRUE),
        paired = metadata$paired,
        threads = metadata$threads
      )

    ccplotpath <- file.path(metadata$outdir, "cross-correlation")
    create_dir(ccplotpath)
    run_cmd(paste("mv *ccplot.png", ccplotpath))
    run_cmd(paste("mv length_results.txt", ccplotpath))
  } else {
    cc <-
      read.table(
        file.path(
          metadata$outdir,
          "cross-correlation",
          "length_results.txt"
        ),
        sep = "\t",
        header = TRUE
      )
  }

  # Add to the metadata file
  bam.files <-
    list.files(file.path(metadata$outdir, "bam"), "*.bam$", full.names = TRUE)
  filtered.bam.files <-
    list.files(file.path(metadata$outdir, "filteredbam"),
               "*.bam$",
               full.names = TRUE)

  df1 <-
    metadata$sampleinfo %>% dplyr::mutate(I = sub("_1$", "", reduce_path(file1)))

  bf <-
    data.frame(
      bam = bam.files,
      I = sub("_1$", "", reduce_path(bam.files)),
      stringsAsFactors = FALSE
    )
  ff <- data.frame(
    filteredbam = filtered.bam.files,
    I = sub("_1$", "", reduce_path(filtered.bam.files)),
    stringsAsFactors = FALSE
  )
  df2 <- dplyr::left_join(df1, bf, by = "I") %>%
    dplyr::left_join(., ff, by = "I") %>%
    dplyr::select(-I)

  sampleinfo <- df2 %>%
    dplyr::mutate(I = sub("_1$", "", reduce_path(input1))) %>%
    dplyr::left_join(., bf, by = "I") %>% dplyr::left_join(., ff, by = "I") %>%
    dplyr::rename(
      bam = bam.x,
      filteredbam = filteredbam.x,
      inputbam = bam.y,
      inputfilteredbam = filteredbam.y
    )

  cc$I <- rownames(cc)
  qc$I <- rownames(qc)
  sampleinfo <-
    dplyr::mutate(sampleinfo, I = sub("_1$", "", reduce_path(filteredbam))) %>%
    dplyr::left_join(., qc, by = "I") %>%
    dplyr::left_join(., cc, by = "I") %>% dplyr::select(-I)

  metadata$sampleinfo <- sampleinfo

  ## STEP 5 (optional): generate track files for genome browser ---------------
  if (bigwig) {
    create_bw(metadata)
    out.bw <- file.path(metadata$outdir, "bw")
    create_dir(out.bw)
    to.move <-
      paste0(reduce_path(metadata$sampleinfo$filteredbam), ".bw")
    mapply(run_cmd, (sprintf("mv %s %s", to.move, out.bw)))
  }

  ## STEP 6: Call peaks on samples --------------------------------------------

  # Use MACS2 to call peaks
  # Note the lax qvalue - this is because we are going to be performing DE analysis
  # using a count-based method.  Using this lax threshold identifies all regions with an
  # appreciable signal

  peak.out <-
    metadata$sampleinfo$sample %>% paste0("_fdr10_peaks.narrowPeak")
  peak.out.path <- file.path(metadata$outdir, "peaks", peak.out)

  if (!all(file.exists(peak.out.path))) {
    give_note("Calling peaks\n")
    if (metadata$paired) {
      if (metadata$input) {
        peaks <-
          run_macs2(
            macs2 = macs2,
            out.dir = file.path(metadata$outdir, "peaks"),
            treatment.files = metadata$sampleinfo$filteredbam,
            qvalue = 0.1,
            control.files = metadata$sampleinfo$inputfilteredbam,
            input.format = "BAMPE",
            species = metadata$annotation,
            sample.names = metadata$sampleinfo$sample,
            call.summits = TRUE,
            threads = metadata$threads,
            broad = broad
          )
      } else {
        peaks <-
          run_macs2(
            macs2 = macs2,
            out.dir = file.path(metadata$outdir, "peaks"),
            treatment.files = metadata$sampleinfo$filteredbam,
            input.format = "BAMPE",
            species = metadata$annotation,
            no.lambda = TRUE,
            no.model = TRUE,
            extsize = metadata$sampleinfo$fragmentlength,
            shift = 0,
            qvalue = 0.1,
            sample.names = metadata$sampleinfo$sample,
            call.summits = TRUE,
            threads = metadata$threads,
            broad = broad
          )
      }
    } else {
      if (metadata$input) {
        peaks <-
          run_macs2(
            macs2 = macs2,
            out.dir = file.path(metadata$outdir, "peaks"),
            treatment.files = metadata$sampleinfo$filteredbam,
            qvalue = 0.2,
            control.files = metadata$sampleinfo$inputfilteredbam,
            input.format = "BAM",
            species = metadata$annotation,
            sample.names = metadata$sampleinfo$sample,
            call.summits = TRUE,
            threads = metadata$threads,
            broad = broad
          )
      } else {
        peaks <-
          run_macs2(
            macs2 = macs2,
            out.dir = file.path(metadata$outdir, "peaks"),
            treatment.files = metadata$sampleinfo$filteredbam,
            qvalue = 0.2,
            input.format = "BAM",
            species = metadata$annotation,
            no.lambda = TRUE,
            no.model = TRUE,
            extsize = metadata$sampleinfo$fragmentlength,
            shift = 0,
            sample.names = metadata$sampleinfo$sample,
            call.summits = TRUE,
            threads = metadata$threads,
            broad = broad
          )
      }
    }
  } else {
    peaks <- lapply(peak.out.path, peak2Granges, broad = broad)
    names(peaks) <- metadata$sampleinfo$sample
  }

  ## STEP 7:  Annotate the peaks and create summary plots ---------------------

  # Merge overlapping peaks and annotate

  if (!file.exists(file.path(metadata$outdir, "peaks", "annotation.rda"))) {
    give_note("Annotating peaks\n")
    annotation <- lapply(peaks, function(x) {
      annotate_peaks(
        x,
        merge = TRUE,
        mergedist = 0,
        tssdist = 1500,
        species = metadata$annotation
      )
    })
    save(annotation,
         file = file.path(metadata$outdir, "peaks", "annotation.rda"))
  } else {
    load(file.path(metadata$outdir, "peaks", "annotation.rda"))
  }

  annotated.peaks <- lapply(annotation, function(x)
    x[[1]])

  null <- lapply(names(annotated.peaks), function(x) {
    df <- Granges2df(annotated.peaks[[x]])
    write.table(
      df,
      file = file.path(
        metadata$outdir,
        "peaks",
        paste0(x, "_annotated_peaks.txt")
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  })

  # Create summary plots
  for (i in seq_along(annotation)) {
    plot_peaksummary(annotation[[i]])
    run_cmd(sprintf(
      "mv peaksDistribution.png %s",
      file.path(
        metadata$outdir,
        "peaks",
        paste0(names(annotation)[[i]], "_peaks_distribution.png")
      )
    ))
  }

  # Create a plot summarising the distance to TSS for all the samples ---------
  dist2TSS <- lapply(annotated.peaks, function(x)
    x$disttotss)
  dist2TSS.cut <-
    lapply(dist2TSS, cut, breaks = c(0, 1e3, 1e4, 1e5, 1e10))
  dist2TSS.table <- sapply(dist2TSS.cut, table)
  dist2TSS.percentage <- apply(dist2TSS.table, 2,
                               function(.ele)
                                 .ele / sum(.ele))
  dist2TSS.percentage <- reshape2::melt(dist2TSS.percentage)
  dist2TSS.percentage$value <- dist2TSS.percentage$value * 100

  dtplot <-
    ggplot2::ggplot(dist2TSS.percentage,
                    ggplot2::aes(x = Var2, y = value, fill = Var1)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("") + ggplot2::ylab("%") +
    ggplot2::labs(title = "Distance to TSS") +
    theme_publication() +
    ggplot2::scale_fill_manual(
      name = "",
      values = c("#FF0000FF", "#FF0000BB",
                 "#FF000077", "#FF000033"),
      labels = c("<1kb", "1-10kb", "10-100kb", ">100kb")
    )

  ggplot2::ggsave(filename = file.path(metadata$outdir, "peaks", "dist2tss.png"),
                  plot = dtplot)


  # For each sample, count reads in peaks to calculate FRIP
  qc.saf <- lapply(names(annotated.peaks), function(x) {
    Granges2saf(annotated.peaks[[x]])
    run_cmd(sprintf("mv annotated.peaks[[x]].saf %s", paste0(x, ".saf")))
    paste0(x, ".saf")
  }) %>%
    unlist()

  count.dir <-
    file.path(metadata$outdir, "peaks", "counts_in_peaks")
  create_dir(count.dir)
  frip <- vector("list", length = length(qc.saf))
  for (i in seq_along(metadata$sampleinfo$filteredbam)) {
    counts <- run_featurecounts(
      featurecounts = featurecounts,
      annotationFile = qc.saf[[i]],
      requireBothEndsMapped = TRUE,
      excludeChimeric = TRUE,
      pairedEnd = metadata$paired,
      countMultiMapping = FALSE,
      multiFeatureReads = FALSE,
      threads = metadata$threads,
      ignoreDup = TRUE,
      outname = paste0(reduce_path(qc.saf[[i]]), ".counts.txt"),
      alignments = metadata$sampleinfo$filteredbam[[i]]
    )
    res <-
      read.table(paste0(reduce_path(qc.saf[[i]]), ".counts.txt.full.summary"),
                 header = TRUE)
    frip[[i]] <- 100 * res[1, 2] / sum(res[, 2])
  }

  frip <- unlist(frip)
  names(frip) <- reduce_path(unlist(qc.saf))
  frip <-
    data.frame(
      sample = names(frip),
      FRIP = unlist(frip),
      stringsAsFactors = FALSE
    )
  sampleinfo <-
    metadata$sampleinfo %>% dplyr::left_join(., frip, by = "sample")
  metadata$sampleinfo <- sampleinfo

  unlink(qc.saf)
  unlink(list.files(pattern = "counts.txt"))

  ## STEP 8: Generate consensus peaks

  # If replicates >2, then peaks must be present in at least half; otherwise present in both
  cond <- metadata$sampleinfo %>%
    split.data.frame(., .$condition, drop = TRUE) %>%
    lapply(., dplyr::pull, sample)
  peaks.by.condition <- lapply(cond, function(x) {
    annotated.peaks[names(annotated.peaks) %in% x]
  })

  nr <- lapply(peaks.by.condition, function(x) {
    if (length(x) <= 2) {
      return(2)
    } else {
      return(ceiling(length(x) / 2))
    }
  })

  con.out <- file.path(metadata$outdir, "consensus")
  create_dir(con.out)

  consensus.peaks <- vector("list", length(peaks.by.condition))
  for (i in seq_along(peaks.by.condition)) {
    consensus.peaks[[i]] <-
      consensus_peaks(peaks.by.condition[[i]],
                      replicates = nr[[i]],
                      plot = TRUE)
    run_cmd(sprintf("mv consensus.png %s",
                    file.path(
                      con.out, paste0(names(peaks.by.condition)[[i]], "_consensus_peaks.png")
                    )))
  }
  names(consensus.peaks) <- names(peaks.by.condition)

  # annotate consensus peaks
  annotation.consensus <- lapply(consensus.peaks, function(x) {
    annotate_peaks(
      x,
      merge = TRUE,
      mergedist = 0,
      tssdist = 1500,
      species = metadata$annotation
    )
  })

  annotated.consensus.peaks <-
    lapply(annotation.consensus, function(x)
      x[[1]])

  null <- lapply(names(annotated.consensus.peaks), function(x) {
    df <- Granges2df(annotated.consensus.peaks[[x]])
    write.table(
      df,
      file = file.path(con.out, paste0(x, "_annotated_peaks.txt")),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  })

  # Create summary plots
  for (i in seq_along(annotation.consensus)) {
    plot_peaksummary(annotation.consensus[[i]])
    run_cmd(sprintf("mv peaksDistribution.png %s",
                    file.path(
                      con.out, paste0(names(annotation.consensus)[[i]],
                                      "_peaks_distribution.png")
                    )))
  }

  # Merge the consensus peaks from the different condition to identify
  # all peaks
  merged.peaks <- GenomicRanges::reduce(consensus.peaks[[1]])
  for (i in seq_along(consensus.peaks)[-1]) {
    merged.peaks <-
      c(merged.peaks, GenomicRanges::reduce(consensus.peaks[[i]]))
  }
  merged.peaks <- GenomicRanges::reduce(merged.peaks)

  # STEP 9: Perform binary differential analysis ------------------------------
  pr <-
    lapply(consensus.peaks, function(x)
      GenomicRanges::findOverlaps(merged.peaks, x))
  dr <-
    data.frame(matrix(
      nrow = length(merged.peaks),
      ncol = length(consensus.peaks),
      dimnames = list(1:length(merged.peaks), 1:length(consensus.peaks))
    ))

  colnames(dr) <- names(consensus.peaks)
  for (i in seq_along(dr)) {
    dr[, i][as.numeric(rownames(dr)) %in% pr[[i]]@from] <- 1
    dr[, i][!(as.numeric(rownames(dr)) %in% pr[[i]]@from)] <- 0
  }

  x <-
    apply(combn(ncol(dr), 2), 2, function(x)
      dr[, x[2]] - dr[, x[1]])
  colnames(x) <- apply(combn(ncol(dr), 2), 2, function(x) {
    paste(rev(names(dr)[x]), collapse = '-')
  })
  GenomicRanges::elementMetadata(merged.peaks) <- x

  diff.binary <-
    file.path(metadata$outdir, "differential", "binary")
  create_dir(diff.binary)

  merged.peaks.df <- Granges2df(merged.peaks)
  write.table(
    merged.peaks.df,
    file = file.path(diff.binary, "binary_results.txt"),
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  # STEP 10: Perform count-based differential analysis ------------------------

  # Create SAF from merged peaks
  saf <- Granges2saf(merged.peaks)

  # Count reads in peak regions
  counts <- run_featurecounts(
    featurecounts = featurecounts,
    annotationFile = "merged.peaks.saf",
    requireBothEndsMapped = TRUE,
    excludeChimeric = TRUE,
    pairedEnd = metadata$paired,
    countMultiMapping = FALSE,
    multiFeatureReads = FALSE,
    threads = metadata$threads,
    ignoreDup = TRUE,
    outname = "counts.txt",
    alignments = metadata$sampleinfo$filteredbam
  )

  metadata$counts <- counts
  unlink(list.files(pattern = "saf"))
  unlink(list.files(pattern = "counts.txt"))

  # Run differential analysis
  tryCatch({
    DE_deseq2 <- deseq2_analysis(metadata = metadata, species = species)
  }, error = function(e)
    cat(crayon::red(
      crayon::bold("\nNo differential peaks found by count-based method\n")
    )))

  if (is.defined(DE_deseq2)) {
    if (nrow(DE_deseq2) > 0)  {
      DE_deseq2 <- lapply(DE_deseq2, function(x) {
        x %>% dplyr::mutate(
          chr = gsub(":.*", "", peak),
          start = sub(".*: *(.*?) *\\-.*", "\\1", peak),
          end = gsub(".*\\-", "", peak)
        ) %>%
          as.data.frame() %>%
          GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      })
    }
    return(list(
      DE_countbased = DE_deseq2,
      DE_binary = merged.peaks.df,
      metadata = metadata
    ))
  } else {
    return(list(DE_binary = merged.peaks.df, metadata = metadata))
  }
}
