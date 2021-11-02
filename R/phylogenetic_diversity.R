#' calculate phylogenetic diversity, e.g. the total branch length of the extant
#' tree
#' @param phy phylo object
#' @param t time point at which to measure phylogenetic diversity
#' @param extinct_tol tolerance to determine if a lineage is extinct at time t.
#' Default is 1/100 * smallest branch length of the tree.
#' @return phylogenetic diversity
#' @export
phylogenetic_diversity <- function(phy,
                                   t = 1e10,
                                   extinct_tol = NULL) {
  crown_age = max(ape::branching.times(phy))[[1]]
  if (is.null(extinct_tol)) {
    extinct_tol <- min(phy$edge.length)/100
  }
  return(calc_phylodiv_cpp(phy, t, crown_age, extinct_tol));
}
