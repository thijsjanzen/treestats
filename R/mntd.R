#' Fast function using C++ to calculate the Mean Nearest Taxon distance
#' @description After calculating all pairwise distances between all tips,
#' this function takes the mean value per tip, and then calculates the average
#' value across all tips.
#' @param phy phylo object or ltable
#' @return Mean Nearest Taxon Distance.
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#' @export
mntd <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_mntd_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_mntd_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
