#' GREAT analysis
#'
#' Wrapper around rGREAT package to perform pathway analysis of genomic regions.
#' Uses the basal-plus-extension rule used as default on the GREAT website.
#'
#' @param regions Granges object of genomic regions to submit for GREAT analysis
#' @param species character string specifying the name of the species. Only
#' \code{'human'}, and \code{'mouse'} are supported at present.  [DEFAULT = human].
#'
#' @export
#' @importFrom GenomeInfoDb mapSeqlevels seqlevels renameSeqlevels
#' @importFrom rtracklayer liftOver
#' @importFrom rGREAT submitGreatJob getEnrichmentTables

run_great <- function(regions = NULL, species = c("human", "mouse")) {

  chain <- NULL
  species <- match.arg(species)

  if (species == "human") {
    utils::data(chain)
    new_style <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(regions), "UCSC")
    reg <- GenomeInfoDb::renameSeqlevels(regions, new_style)
    reg <- unlist(rtracklayer::liftOver(reg, chain))
    build <- "hg19"
  }

  if (species == "mouse") {
    new_style <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(regions), "UCSC")
    reg <- GenomeInfoDb::renameSeqlevels(regions, new_style)
    build <- "mm10"
  }

  job <- rGREAT::submitGreatJob(gr = reg, species = build, request_interval = 30)
  tb <- rGREAT::getEnrichmentTables(job, category = c("GO", "Pathway Data"))
  return(tb)
}
