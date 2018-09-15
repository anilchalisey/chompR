#' File conversions
#'
#' @description Function to convert a narrowPeak or broadPeak MACS2 file into a
#'   GRanges object
#'
#' @param peakfile Path to MACS2 peak file in narrowPeak or broadPeak format
#' @param broad Logical specifying whether broad or narrowpeak
#'
#' @return
#' GRanges object of peak regions
#'
#' @export
#' @importFrom rtracklayer import
#'
#' @examples
#' \dontrun{
#' peaks <- peak2Granges("peaks/sample.narrowPeak")
#' }

peak2Granges <- function(peakfile, broad = FALSE) {
  extraCols.narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  extraCols.broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                           qValue = "numeric")
  broad <- ifelse(grepl(".broadPeak", peakfile), TRUE, FALSE)

  if (isTRUE(broad)) {
    grg <- rtracklayer::import(peakfile, format = "BED",
                               extraCols = extraCols.broadPeak)
  } else {
    grg <- rtracklayer::import(peakfile, format = "BED",
                               extraCols = extraCols.narrowPeak)
  }
  return(grg)
}


#' @description Function to convert a GRanges object to a dataframe
#'
#' @rdname peak2Granges
#'
#' @param gr GRanges object to convert
#' @param ignore.strand Logical indicating whether strand information should be
#'   included
#'
#' @return
#' A dataframe created from the GRanges object
#'
#' @export
#'
#' @importFrom GenomicRanges seqnames start end strand elementMetadata
#' @importFrom S4Vectors DataFrame

Granges2df <- function(gr, ignore.strand = FALSE) {
  df.columns <- list(row.names = names(gr),
                     chr = GenomicRanges::seqnames(gr),
                     start = GenomicRanges::start(gr),
                     end = GenomicRanges::end(gr))
  if (!ignore.strand) {
    df.columns <- c(df.columns, list(strand = GenomicRanges::strand(gr)))
  }
  df.columns <- c(df.columns, as.list(GenomicRanges::elementMetadata(gr)))
  do.call(S4Vectors::DataFrame, c(df.columns, check.names = FALSE))
}

#' @description Function to convert a GRanges object to a SAF
#'
#' @rdname peak2Granges
#'
#' @return
#' A SAF file created from the GRanges object
#'
#' @export
#' @importFrom S4Vectors DataFrame

Granges2saf <- function(gr) {
  GeneID <- paste(GenomicRanges::seqnames(gr), GenomicRanges::start(gr), sep = ":")
  GeneID <- paste(GeneID, GenomicRanges::end(gr), sep = "-")
  saf.columns <- list(GeneID = GeneID,
                      Chr = GenomicRanges::seqnames(gr),
                      Start = GenomicRanges::start(gr),
                      End = GenomicRanges::end(gr),
                      Strand = GenomicRanges::strand(gr))
  saf <- do.call(S4Vectors::DataFrame, saf.columns)
  saf.name <- paste0(deparse(substitute(gr)), ".saf")
  write.table(saf, file = saf.name, row.names = FALSE, col.names = TRUE,
              quote = FALSE, sep = "\t")
  return(saf)
}
