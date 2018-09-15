#' Annotation of peaks
#'
#' This function categorises peaks as either promoter, genebody or intergenic
#' based on the distance from a known transcription start site (TSS) and gene
#' overlap.  By default, the distance is set to 1500 bp but can be changed with
#' the 'tssdist' argument.  The function also allows regions that are within a
#' certain distance of each other to be merged to form a larger region.  The
#' merge distance is set as default to 0 (i.e. they have to be touching) but may be
#' changed with the 'mergedist' argument, and will only merge peaks within the
#' same region (i.e. distal peaks will not be merged with proximal peaks).
#'
#' @param regions GRanges object of peaks to annotate
#' @param merge Boolean indicating whether to merge peaks or not
#' @param mergedist Maximum distance apart for peaks to be merged
#' @param tssdist Distance from TSS to classify as a proximal or distal peak
#' @param species character string specifying the name of the species. Only
#' \code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].
#'
#' @return
#' GRanges object of annotated peaks
#'
#' @export
#'
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom TxDb.Mmusculus.UCSC.mm10.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom BiocInstaller biocLite
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom GenomicRanges start end reduce sort elementMetadata
#' @importFrom ChIPpeakAnno assignChromosomeRegion annoPeaks annoGR
#' @importFrom dplyr left_join select
#' @importFrom AnnotationDbi toTable
#' @importFrom BiocInstaller biocLite

annotate_peaks <- function(regions, merge = TRUE,
                            mergedist = 0,
                            tssdist = 1500,
                            species = c("human", "mouse")) {

  # function to annotate peaks as proximal or distal to tss

  disttotss <- region <- feature <- symbol <- NULL

  species <- match.arg(species)

  if (species == "human") {
    org <- "Homo sapiens"
    build <- "hg38"
    slevels <- paste0("chr", c(1:22, "X", "Y"))
  }

  if (species == "mouse") {
    org <- "Mus musculus"
    build <- "mm10"
    slevels <- paste0("chr", c(1:19, "X", "Y"))
  }

  txdb <- paste0("TxDb.", toupper(substr(org, 1, 1)),
                 gsub("^.* ", "", org), ".UCSC.", build,
                 ".knownGene")

  if (!requireNamespace(txdb, quietly = TRUE)) BiocInstaller::biocLite(txdb)
  txdb <- eval(parse(text = paste0(txdb, "::", txdb)))
  GenomeInfoDb::seqlevels(txdb) <- slevels
  tx <- ChIPpeakAnno::annoGR(txdb)
  annotated_regions <- tss_granges(reg = regions, tssdist = tssdist, tx = tx)

  if (merge) {
    merged_regions <- annotated_regions
    GenomicRanges::start(merged_regions) <-
      GenomicRanges::start(merged_regions) - round(mergedist/2)
    GenomicRanges::end(merged_regions) <-
      GenomicRanges::end(merged_regions) + round(mergedist/2)
    prox_peaks <- merged_regions[merged_regions$region == "proximal"]
    dist_peaks <- merged_regions[merged_regions$region == "distal"]
    prox_peaks <- GenomicRanges::reduce(prox_peaks)
    dist_peaks <- GenomicRanges::reduce(dist_peaks)
    GenomicRanges::start(prox_peaks) <-
      GenomicRanges::start(prox_peaks) + round(mergedist/2)
    GenomicRanges::end(dist_peaks) <-
      GenomicRanges::end(dist_peaks) - round(mergedist/2)
    merged_regions <- GenomicRanges::sort(c(prox_peaks, dist_peaks))
    annotated_regions <- tss_granges(merged_regions, tssdist, tx)
  }

  positions <- ChIPpeakAnno::assignChromosomeRegion(merged_regions, nucleotideLevel = FALSE,
                                                    precedence = c("Promoters", "immediateDownstream",
                                                                   "fiveUTRs", "threeUTRs", "Exons", "Introns"),
                                                    TxDb = txdb, proximal.promoter.cutoff = 1000,
                                                    immediate.downstream.cutoff = tssdist)

  anno_regions <- ChIPpeakAnno::annoPeaks(annotated_regions, annoData = tx, bindingType = "startSite",
                                          select = "bestOne", bindingRegion = c(-125000000, 125000000))

  orgdb <- paste0("org.", toupper(substr(org, 1, 1)),
                  tolower(substr(gsub("^.* ", "", org), 1, 1)),
                  ".eg.db")
  symb <- paste0(gsub(".db", "", orgdb), "SYMBOL")
  orgdb <- eval(parse(text = paste0(orgdb, "::", symb)))
  e2s <- AnnotationDbi::toTable(orgdb)
  colnames(e2s) <- c("feature", "symbol")
  em <- as.data.frame(GenomicRanges::elementMetadata(anno_regions))
  em <- dplyr::left_join(em, e2s) %>% dplyr::select(disttotss, region, entrezid = feature, symbol)
  GenomicRanges::elementMetadata(anno_regions) <- em
  list(peaks = anno_regions, peak.regions = positions)
}

#' Classify regions as either proximal or distal to TSS
#'
#' Helper function for annotate peaks
#'
#' @param reg GRanges object of peaks to annotate
#' @param tssdist Distance from TSS to classify as a proximal or distal peak
#' @param edb ensemblDB object containing Ensembl annotation for relevant genome
#'
#' @noRd

tss_granges <- function(reg, tssdist, tx) {
  new_style <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(reg), "UCSC")
  reg <- GenomeInfoDb::renameSeqlevels(reg, new_style)
  disttotss <- GenomicRanges::distanceToNearest(
    reg, GenomicRanges::reduce(GenomicRanges::resize(tx, 1, fix = "start")))
  GenomicRanges::mcols(reg)$disttotss <-
    GenomicRanges::mcols(disttotss)$distance
  GenomicRanges::mcols(reg)$transcript_id <-
    GenomicRanges::mcols(tx[S4Vectors::subjectHits(disttotss), ])$transcript_id
  GenomicRanges::mcols(reg)$gene_id <-
    GenomicRanges::mcols(tx[S4Vectors::subjectHits(disttotss), ])$gene_id
  GenomicRanges::mcols(reg)$gene_name <-
    GenomicRanges::mcols(tx[S4Vectors::subjectHits(disttotss), ])$gene_name
  GenomicRanges::mcols(reg)$region <-
    ifelse(GenomicRanges::mcols(reg)$disttotss >= tssdist,
           "distal", "proximal")
  return(reg)
}
