#' Fast function using C++ to calculate the sackin index of (im)balance.
#' @description The Sackin index is calculated as the sum of ancestors for each
#' of the tips. Higher values indicate higher imbalance. Two normalizations
#' are available, where a correction is made for tree size, under either a
#' Yule expectation, or a pda expectation.
#' @param phy phylogeny or ltable
#' @return Sackin index
#' @references M. J. Sackin (1972). "Good" and "Bad" Phenograms.
#' Systematic Biology. 21:225-226.
#' @export
var_leaf_depth <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_sackin_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_var_leaf_depth_cpp(as.vector(t(phy$edge))))
  }
}
