#' Sackin index of (im)balance.
#' @description The Sackin index is calculated as the sum of ancestors for each
#' of the tips. Higher values indicate higher imbalance. Two normalizations
#' are available, where a correction is made for tree size, under either a
#' Yule expectation, or a pda expectation.
#' @param phy phylogeny or ltable
#' @param normalization normalization, either 'none' (default), "yule" or "pda".
#' @return Sackin index
#' @references M. J. Sackin (1972). "Good" and "Bad" Phenograms.
#' Systematic Biology. 21:225-226.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#' unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#' sackin(balanced_tree)
#' sackin(unbalanced_tree) # should be much higher
sackin <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(phy, "matrix")) {
    return(calc_sackin_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_sackin_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}
