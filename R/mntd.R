#' Mean Nearest Taxon distance
#' @description Per tip, evaluates the shortest distance to another tip,
#' then takes the average across all tips.
#' @param phy phylo object or ltable
#' @return Mean Nearest Taxon Distance.
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#' @export
mntd <- function(phy) {

  # only defined for extant trees!
  if (inherits(phy, "matrix")) {
    if (sum(phy[, 4] != -1)) {
      stop("can only calculate mntd statistic for ultrametric tree")
    }
    return(calc_mntd_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    if (!ape::is.ultrametric(phy, option = 2)) {
      stop("can only calculate mntd statistic for ultrametric tree")
    }
    return(calc_mntd_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
