#' Root imbalance
#' @description Measures the distribution of tips over the two crown lineages,
#' e.g. n1 / (n1 + n2), where n1 is the number of tips connected to crown
#' lineage 1 and n2 is the number of tips connected to crown lineage 2. We
#' always take n1 > n2, thus root imbalance is always in [0.5, 1].
#' @param phy phylo object or ltable
#' @return Root imbalance
#' @references  Guyer, Craig, and Joseph B. Slowinski. "Adaptive radiation and
#' the topology of large phylogenies." Evolution 47.1 (1993): 253-263.
#' @export
root_imbalance <- function(phy) {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(phy, "matrix")) {
    return(calc_root_imbalance_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_root_imbalance_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}
