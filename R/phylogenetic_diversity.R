#' calculate phylogenetic diversity, e.g. the total branch length of the extant
#' tree
#' @param phy phylo object
#' @param t time point at which to measure phylogenetic diversity
#' @return phylogenetic diversity
#' @export
phylogenetic_diversity <- function(phy,
                                   t = 1e10) {
  has_extinct <- length(geiger::is.extinct(phy)) > 0
  crown_age = max(ape::branching.times(phy))[[1]]
  return(calc_phylodiv_cpp(phy, t, has_extinct, crown_age));
}
