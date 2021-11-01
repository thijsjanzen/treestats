#' calculate colless index, using apTreeshape
#' @param phy phylo object
#' @param norm A character string equals to NULL (default) for no
#' normalization or one of "pda" or "yule".
#' @return colless index
#' @export
colless <- function(phy,
                    norm = NULL) {
  if (class(phy) != "phylo") {
    stop("input has to be phylo object")
  }
  return(apTreeshape::colless(apTreeshape::as.treeshape(phy), norm))
}
