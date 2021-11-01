#' calculate phylogenetic diversity, e.g. the total branch length of the extant
#' tree
#' @param phy phylo object
#' @param t time point at which to measure phylogenetic diversity
#' @return phylogenetic diversity
#' @export
phylogenetic_diversity <- function(phy,
                                   t = 1e10) {
  return(calc_phylodiv_cpp(phy, t));
}
