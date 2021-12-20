#' calculate the nLTT, using C++. Only use if you are very certain
#' about the input data, and are certain that performing nLTT is valid (e.g.
#' your tree is ultrametric etc). If you are less certain, use the nLTT function
#' from the nLTT package. The treestats implementation is approximately N/20
#' times faster, where N is the number of extant tips in the largest tree.
#' @param phy phylo or multiPhylo object
#' @param ref_tree reference tree to compare with
#' @return number of lineages
#' @export
nLTT <- function(phy, # nolint
                 ref_tree) {

  brts1 <- treestats::branching_times(phy)
  brts1 <- c(-1.0 * rev(sort(brts1)), 0)
  brts2 <- treestats::branching_times(ref_tree)
  brts2 <- c(-1.0 * rev(sort(brts2)), 0)
  return(calc_nltt_cpp(brts1, brts2))
}

#' calculate the nLTT, using a reference 'empty' tree with only two lineages.
#' Only use if you are very certain about the input data.
#' If you are less certain, use the nLTT function from
#' the nLTT package.
#' @param phy phylo or multiPhylo object
#' @return number of lineages
#' @export
nLTT_base <- function(phy) {  # nolint
  empty_tree <- ape::read.tree(text = "(1:4,2:4):0;")
  return(nLTT(phy, empty_tree))
}
