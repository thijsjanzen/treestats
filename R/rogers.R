#' Rogers J index of (im)balance.
#' @description The Rogers index is calculated as the total number of internal
#' nodes that are unbalanced, e.g. for which both daughter nodes lead to a
#' different number of extant tips. in other words, the number of nodes where
#' L != R (where L(R) is the number of extant tips of the Left (Right) daughter
#' node).
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", in which case the resulting statistic
#' is divided by the number of tips - 2 (e.g. the maximum value of the rogers
#' index for a tree).
#' @return Rogers index
#' @references  J. S. Rogers. Central Moments and Probability Distributions of
#' Three Measures of Phylogenetic Tree Imbalance. Systematic Biology,
#' 45(1):99-110, 1996. doi: 10.1093/sysbio/45.1.99.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#' unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#' rogers(balanced_tree)
#' rogers(unbalanced_tree) # should be higher
rogers <- function(phy, normalization = "none") {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    rogers_stat <- calc_rogers_ltable_cpp(phy)
    if (normalization == "tips" || normalization == TRUE) {
      rogers_stat <- rogers_stat / (length(phy[, 1]) - 2)
    }
    return(rogers_stat)
  }
  if (inherits(phy, "phylo")) {
    rogers_stat <- calc_rogers_cpp(as.vector(t(phy$edge)))
    if (normalization == "tips" || normalization == TRUE) {
      rogers_stat <- rogers_stat / (length(phy$tip.label) - 2)
    }
    return(rogers_stat)
  }
  stop("input object has to be phylo or ltable")
}
