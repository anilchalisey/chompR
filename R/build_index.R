#' Build index for HISAT2
#'
#' \code{build_index} for mapping reads using HISAT2. This actually downloads the
#' pre-compiled index available on the HISAT2 website:
#' \url{ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data}
#'
#' @inheritParams run_chip
#'
#' @export

build_index <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (species == "human") {
    run_cmd("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz")
    idx <- "grch38.tar.gz"
  }
  if (species == "mouse") {
    run_cmd("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38.tar.gz")
    idx <- "grcm38.tar.gz"
  }

  run_cmd(sprintf("tar -xvzf %s", idx))
  run_cmd(sprintf("rm -rf %s", idx))
  file.path(reduce_path(idx), "genome")
}
