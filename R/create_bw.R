#' Create bigwig files
#'
#' @param metadata data.frame object created by run_chip or run_atac funcitons
#'
#' @export
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom TxDb.Mmusculus.UCSC.mm10.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene
#' @importFrom ATACseqQC readBamFile
#' @importFrom GenomeInfoDb seqlevelsStyle mapSeqlevels seqlevels renameSeqlevels
#' @importFrom GenomicAlignments coverage seqinfo
#' @importFrom rtracklayer export.bw

create_bw <- function(metadata) {
  species <- metadata$annotation
  if (species == "human") txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  if (species == "mouse") txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  cvgs <- mapply(function(x, y) {
    gal <- ATACseqQC::readBamFile(x)
    cvg <- GenomicAlignments::coverage(x, width = as.numeric(y))
    new_style <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(cvg), "UCSC")
    reg <- GenomeInfoDb::renameSeqlevels(cvg, new_style)
    GenomicAlignments::seqinfo(cvg) <- GenomicAlignments::seqinfo(txdb)
    rtracklayer::export.bw(cvg, con = paste0(reduce_path(x), ".bw"))
    cvg
  }, metadata$sampleinfo$filteredbam, metadata$sampleinfo$fragmentlength
  )
}
