#' Calculate the nLTT, using C++.
#' @description The nLTT statistic calculates the sum of
#' absolute differences in the number of lineages over time, where both the
#' number of lineages and the time are normalized. The number of lineages is
#' normalized by the number of extant tips, whereas the time is normalized by
#' the crown age. The nLTT can only be calculated for reconstructed trees.
#' Only use the treestats version if you are very certain
#' about the input data, and are certain that performing nLTT is valid (e.g.
#' your tree is ultrametric etc). If you are less certain, use the nLTT function
#' from the nLTT package.
#' @param phy phylo object or ltable
#' @param ref_tree reference tree to compare with (should be same type as phy)
#' @references Janzen, T., Höhna, S. and Etienne, R.S. (2015), Approximate
#' Bayesian Computation of diversification rates from molecular phylogenies:
#' introducing a new efficient summary statistic, the nLTT. Methods Ecol Evol,
#' 6: 566-575. https://doi.org/10.1111/2041-210X.12350
#' @return number of lineages
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' reference_tree <- ape::rphylo(n = 10, birth = 0.2, death = 0)
#' nLTT(simulated_tree, reference_tree)
#' nLTT(simulated_tree, simulated_tree) # should be zero.
nLTT <- function(phy, # nolint
                 ref_tree) {

  if (inherits(phy, "phylo")) {
    return(calc_nltt_cpp(phy, ref_tree))
  }

  if (inherits(phy, "matrix")) {
    if (inherits(ref_tree, "phylo")) {
       ref_tree <- treestats::phylo_to_l(ref_tree)
    }
    return(calc_nltt_ltable_cpp(phy, ref_tree))
  }

  stop("input needs to be phylo or ltable object")
}

#' Calculates the nLTT statistic using a reference 'empty' tree with only
#' two lineages.
#' @description The base nLTT statistic can be used as a semi stand-alone
#' statistic for phylogenetic trees. However, please note that although this
#' provides a nice way of checking the power of the nLTT statistic without
#' directly comparing two trees, the nLTT_base statistic is not a substitute
#' for directly comparing two phylogenetic trees. E.g. one would perhaps
#' naively assume that \eqn{nLTT(A, B) = |nLTT(A, base) - nLTT(B, base)}.
#' Indeed, in some cases this may hold true (when, for instance, all normalized
#' lineages of A are less than all normalized lineages of B), but once the
#' nLTT curve of A intersects the nLTT curve of B, this no longer applies.
#' @param phy phylo object
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' nLTT_base(simulated_tree)
#' @return number of lineages
#' @export
nLTT_base <- function(phy) {  # nolint
  empty_tree <- ape::read.tree(text = "(1:4,2:4):0;")
  return(nLTT(phy, empty_tree))
}
