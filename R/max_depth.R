#' Fast function using C++ to calculate maximum depth metric
#' @description The maximum depth metric, measures the maximal path (in edges),
#' between the tips and the root.
#' @param phy phylo object or ltable
#' @return Rogers index
#' @references  C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' @export
max_depth <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_max_depth_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_max_depth_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}
