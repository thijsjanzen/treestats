#' Fast function using C++ to calculate the Blum index of (im)balance.
#' @description The Blum index of imbalance (also known as the s-shape
#' statistic) calculates the sum of \eqn{log(N-1)} over all internal nodes,
#' where N represents the total number of extant tips connected to that node.
#' An alternative implementation can be found in the Castor R package.
#' @param phy phylogeny or ltable
#' @param normalization because the Blum index sums over all nodes,
#' the resulting statistic tends to be correlated with the number of extant
#' tips. Normalization can be performed by dividing by the number of extant
#' tips.
#' @return Blum index of imbalance
#' @references M. G. B. Blum and O. Francois (2006). Which random processes
#' describe the Tree of Life? A large-scale study of phylogenetic tree
#' imbalance. Systematic Biology. 55:685-691.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' if (requireNamespace("nodeSub")) {
#'   balanced_tree <- nodeSub::create_balanced_tree(brts)
#'   unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#'   blum(balanced_tree)
#'   blum(unbalanced_tree) # should be higher
#' }
blum <- function(phy,
                 normalization = FALSE) {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    return(calc_blum_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_blum_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}
