#' Fast function using C++ to calculate the Colless index of (im)balance.
#' @description The Colless index is calculated as the sum of
#' \eqn{abs(L - R)} over all nodes, where L (or R) is the number of extant tips
#' associated with the L (or R) daughter branch at that node.  Higher values
#' indicate higher imbalance. Two normalizations are available,
#' where a correction is made for tree size, under either a yule expectation,
#' or a pda expectation.
#' @param phy phylo object or ltable
#' @param normalization A character string equals to "none" (default) for no
#' normalization or one of "pda" or "yule".
#' @return colless index
#' @references  Colless D H. 1982. Review of: Phylogenetics: The Theory and
#' Practice of Phylogenetic Systematics. Systematic Zoology 31:100-104.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' if (requireNamespace("nodeSub")) {
#'   balanced_tree <- nodeSub::create_balanced_tree(brts)
#'   unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#'   colless(balanced_tree)
#'   colless(unbalanced_tree) # should be higher
#' }
colless <- function(phy,
                    normalization = "none") {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    return(calc_colless_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_colless_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}

#' Fast function using C++ to calculate the equal weights Colless index of
#' (im)balance.
#' @description The equal weights Colless index is calculated as the sum of
#' \eqn{abs(L - R) / (L + R - 2)} over all nodes where L + R > 2,
#' where L (or R) is the number of extant tips associated with the L (or R)
#' daughter branch at that node.  Maximal imbalance is associated with a value
#' of 1.0. The ew_colless index is not sensitive to tree size.
#' @param phy phylo object or ltable
#' @return colless index
#' @references  A. O. Mooers and S. B. Heard. Inferring Evolutionary Process
#' from Phylogenetic Tree Shape. The Quarterly Review of Biology, 72(1), 1997.
#' doi: 10.1086/419657.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' if (requireNamespace("nodeSub")) {
#'   balanced_tree <- nodeSub::create_balanced_tree(brts)
#'   unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#'   ew_colless(balanced_tree)
#'   ew_colless(unbalanced_tree) # should be higher
#' }
ew_colless <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_eWcolless_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_eWcolless_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}
