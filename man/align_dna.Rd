% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align_dna.R
\name{align_dna}
\alias{align_dna}
\title{Alignment of DNA-seq reads}
\usage{
align_dna(threads = 1, output.dir = ".", hisat2 = "hisat2",
  samtools = "samtools", sambamba = "sambamba", species = c("human",
  "mouse"), idx = NULL, reads1 = NULL, reads2 = NULL, fastq = TRUE,
  fasta = FALSE, softClipPenalty = NULL, noSoftClip = FALSE,
  tmo = FALSE, secondary = FALSE, maxAlign = NULL, nomixed = FALSE,
  nodiscordant = FALSE, rgid = NULL, quiet = FALSE,
  non_deterministic = FALSE, maxInsert = NULL, memory = "1G",
  remove.mitochondrial = "MT", remove.duplicates = TRUE,
  hash_table = 262144, overflow_size = 2e+05, io_buffer = 128)
}
\arguments{
\item{threads}{an integer value indicating the number of parallel threads to
be used by FastQC. [DEFAULT = maximum number of available threads - 1].}

\item{output.dir}{character string specifying the directory to which results
will be saved.  If the directory does not exist, it will be created.}

\item{hisat2}{a character string specifying the path to the hisat2 executable.
[DEFAULT = "hisat2"].}

\item{samtools}{a character string specifying the path to the samtools executable.
[DEFAULT = "samtools"].}

\item{sambamba}{a character string specifying the path to the sambamba executable.
[DEFAULT = "sambamba"].}

\item{species}{character string specifying the name of the species. Only
\code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].}

\item{idx}{character vector specifying the basename of the index for the reference genome.
The basename is the name of any of the index files up to but not including the final
.1.ht2, etc. If \code{NULL} then the index for the relevant species (human or mouse) will be
created using the \code{build_index()} function.}

\item{reads1}{Character vector of mate1 reads.  If specified, then reads.dir
must be NULL.}

\item{reads2}{Character vector of mate2 reads.  If specified, then reads.dir
must be NULL.  Must be the same length as mate1.  If single-end sequencing,
then should be left as NULL.}

\item{fastq}{Logical indicating if reads are FASTQ files.}

\item{fasta}{Logical indicating if reads are FASTA files.}

\item{softClipPenalty}{Sets the maximum (MX) and minimum (MN) penalties for
soft-clipping per base, both integers.  Must be given in the format "MX,MN".}

\item{noSoftClip}{Logical indicating whether to disallow soft-clipping.}

\item{tmo}{Logical indicating whether to report only those reads aligning to
known transcripts.}

\item{secondary}{Logical indicating whether to report secondary alignments.}

\item{maxAlign}{Integer indicating the maximum number of distinct primary
alignments to search for each read.}

\item{nomixed}{By default, when hisat2 cannot find a concordant or discordant
alignment for a pair, it then tries to find alignments for the individual
mates. If TRUE, this option disables that behavior.}

\item{nodiscordant}{By default, hisat2 looks for discordant alignments if it
cannot find any concordant alignments.  If true, this option disables that
behavior.}

\item{rgid}{Character string, to which the read group ID is set.}

\item{quiet}{If TRUE, print nothing except alignments and serious errors.}

\item{non_deterministic}{When set to TRUE, HISAT2 re-initializes its
pseudo-random generator for each read using the current time.}

\item{maxInsert}{The maximum fragment length for valid paired-end alignments.
This option is valid only with noSplice = TRUE.}

\item{memory}{String specifying maximum memory per thread; suffix K/M/G
recognized.}

\item{remove.mitochondrial}{Character string.  If set, this will remove reads
mapping to the mitochondrial genome.  The string should match the reference
name for the mitochindrial genome in the alignment file.  Examples include
"ChrM", "M" and "MT".}

\item{remove.duplicates}{If TRUE, duplicate reads will be removed.}

\item{hash_table}{Size of hash table for finding read pairs (default is
262144 reads); will be rounded down to the nearest power of two.  For best
performance should be > (average coverage) * (insert size).}

\item{overflow_size}{Size of the overflow list where reads, thrown out of
the hash table, get a second chance to meet their pairs (default is
200000 reads); increasing the size reduces the number of temporary files
created.}

\item{io_buffer}{Controls sizes of the two buffers (in MB) used for reading
and writing BAM during the second pass (default is 128).}
}
\value{
Raw and filtered BAM files
}
\description{
Alignment of DNA-seq reads
}
