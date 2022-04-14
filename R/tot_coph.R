#' Fast function using C++ to calculate the total cophenetic index.
#' @description The total cophenetic index is the sum of the depth of the last
#' common ancestor of all pairs of leaves.
#' @param phy phylo object or ltable
#' @return Total cophenetic index
#' @references  A. Mir, F. Rosselló, and L. Rotger. A new balance index for
#' phylogenetic trees. Mathematical Bio-sciences, 241(1):125-136, 2013.
#' doi: 10.1016/j.mbs.2012.10.005.
#' @export
tot_coph <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_tot_coph_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_tot_coph_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}
