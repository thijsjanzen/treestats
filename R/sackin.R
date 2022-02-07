#' Fast function using C++ to calculate the sackin index of (im)balance. The
#' sackin index is calculated as the sum of ancestors for each of the tips.
#' Higher values indicate higher imbalance. Two normalizations are available,
#' where a correction is made for tree size, under either a yule expectation,
#' or a pda expectation. The sackin index is only available in treestats
#' for extant, ultrametric, strictly bifurcating, trees
#' @param phy phylogeny or ltable
#' @param normalization normalization, either 'none' (default), "yule" or "pda".
#' @return sackin index
#' @references M. J. Sackin (1972). "Good" and "Bad" Phenograms.
#' Systematic Biology. 21:225-226.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' balanced_tree <- nodeSub::create_balanced_tree(brts)
#' unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#' sackin(balanced_tree)
#' sackin(unbalanced_tree) # should be much higher
sackin <- function(phy, normalization = "none") {
  if (inherits(phy, "phylo")) {
    if (!ape::is.ultrametric(phy)) {
      stop("can only calculate beta statistic for ultrametric tree")
    }

    sackin_index <- calc_sackin_cpp(phy, normalization)
    return(sackin_index)
  }

  if (inherits(phy, "matrix")) {
    sackin_index <- calc_sackin_ltable_cpp(phy, normalization)
    return(sackin_index)
  }

  stop("input object has to be phylo or ltable")
}
