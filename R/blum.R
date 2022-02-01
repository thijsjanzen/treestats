#' Fast function using C++ to calculate the Blum index of (im)balance. The Blum
#' index of imbalance calculates the sum of \eqn{log(N-1)} over all internal
#' nodes, where N represents the total number of extant tips connected to that
#' node. The Blum statistic thus is only available for ultrametric,
#' reconstructed trees. An alternative implementation can be found in the
#' Castor R package.
#' @param phy phylogeny or ltable
#' @return Blum index of imbalance
#' @references M. G. B. Blum and O. Francois (2006). Which random processes
#' describe the Tree of Life? A large-scale study of phylogenetic tree
#' imbalance. Systematic Biology. 55:685-691.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' balanced_tree <- nodeSub::create_balanced_tree(brts)
#' unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#' blum(balanced_tree)
#' blum(unbalanced_tree) # should be higher
blum <- function(phy) {
  if (!ape::is.ultrametric(phy)) {
    stop("can only calculate Blum statistic for ultrametric tree")
  }
  blum_index <- apply_function_phy(phy, calc_blum_cpp)
  return(blum_index)
}
