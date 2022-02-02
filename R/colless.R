#' Fast function using C++ to calculate the colless index of (im)balance. The
#' colless index is calculated as the sum of \eqn{abs(L - R)} over all nodes,
#' where L (or R) is the number of extant tips associated with the L (or R)
#' daughter branch at that node
#' Higher values indicate higher imbalance. Two normalizations are available,
#' where a correction is made for tree size, under either a yule expectation,
#' or a pda expectation. The sackin index is only available in treestats
#' for extant, ultrametric, strictly bifurcating, trees. For trees including
#' extinct species, we advise to use the slower version in the package
#' apTreeshape or Castor.
#' @param phy phylo object
#' @param normalization A character string equals to NULL (default) for no
#' normalization or one of "pda" or "yule".
#' @return colless index
#' @references  Colless D H. 1982. Review of: Phylogenetics: The Theory and
#' Practice of Phylogenetic Systematics. Systematic Zoology 31:100-104.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' balanced_tree <- nodeSub::create_balanced_tree(brts)
#' unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#' colless(balanced_tree)
#' colless(unbalanced_tree) # should be higher
colless <- function(phy,
                    normalization = "none") {
  if (!ape::is.ultrametric(phy)) {
    stop("can only calculate colless statistic for ultrametric tree")
  }

  colless_index <- apply_function_phy(phy, calc_colless_cpp, normalization)
  return(colless_index)
}
