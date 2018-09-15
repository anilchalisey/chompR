#' Shift ATAC-seq reads
#'
#' Tn5 transposase ahs been shown to bind as a dimer and inserts two adaptors into
#' accessible DNA locations separated by 9bp.  This function offsets reads on the
#' positive strand by +4bp and those on the negative strand by -5bp to account for
#' this for downstream analysis.
#'
#' @param bam.file character string specifying path to bam file
#' @param species character string specifying the name of the species. Only
#' \code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].
#' @param paired logical, if \code{TRUE} then paired end reads
#'
#' @export
#' @importFrom ATACseqQC readBamFile shiftGAlignmentsList
#' @importFrom rtracklayer export

shift_bam <- function(bam.file, species = c("human", "mouse"), paired = TRUE) {
  tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
  out <- "shifted"
  create_dir(out)
  gal <- ATACseqQC::readBamFile(bam.file, tag = tags, asMates = paired)
  gal1 <- ATACseqQC::shiftGAlignmentsList(gal)
  rtracklayer::export(gal1, file.path(out, paste0(reduce_path(bam.file), ".bam")))
}
