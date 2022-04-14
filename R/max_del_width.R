#' Fast function using C++ to calculate the maximum difference of widths of a
#' phylogenetic tree
#' @description Calculates the maximum difference of widths of a phylogenetic
#' tree. First, the widths are calculated by collecting the depth of each node
#' and tip across the entire tree, where the depth represents the distance
#' (in nodes) to the root. Then, the width represents the number of occurrences
#' of each possible depth. Then, we take the difference between each consecutive
#' width, starting with the first width. The maximum difference is then
#' returned.
#' @param phy phylogeny or ltable
#' @return maximum difference of widths
#' @references C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018..
#' @export
max_del_width <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_max_del_width_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_max_del_width_cpp(as.vector(t(phy$edge))))
  }
}
