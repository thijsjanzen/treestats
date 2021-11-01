#' calculate the nLTT, using C++. Only use if you are very certain
#' about the input data. If you are less certain, use the nLTT function from
#' the nLTT package.
#' @param phy phylo or multiPhylo object
#' @param ref_tree reference tree to compare with
#' @return number of lineages
#' @export
nLTT <- function(phy, # nolint
                 ref_tree) {

  brts1 <- ape::branching.times(phy)
  brts1 <- c(-1.0 * rev(sort(brts1)), 0)
  brts2 <- ape::branching.times(ref_tree)
  brts2 <- c(-1.0 * rev(sort(brts2)), 0)
  return(calc_nltt_cpp(brts1, brts2))
}

#' calculate the nLTT, using a reference 'empty' tree with only two lineages.
#' Only use if you are very certain
#' about the input data. If you are less certain, use the nLTT function from
#' the nLTT package.
#' @param phy phylo or multiPhylo object
#' @return number of lineages
#' @export
nLTT_base <- function(phy) {  # nolint
  empty_tree <- ape::read.tree(text = "(1:4,2:4):0;")
  return(nLTT(phy, empty_tree))
}
