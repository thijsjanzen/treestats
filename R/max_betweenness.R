#' Fast function using C++ to calculate maximum betweenness centrality.
#' @description Betweenness centrality associates with each node v, the two
#' nodes u, w, for which the shortest path between u and w runs through v, if
#' the tree were re-rooted at node v. Then, we report the node with maximum
#' betweenness centrality.
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", if tips is chosen, the obtained
#' maximum betweenness is normalized by the total amount of node pair
#' combinations considered, e.g. (n-2)*(n-1), where n is the number of tips.
#' @return Maximum Betweenness
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
max_betweenness <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    betweenness_stat <- calc_max_betweenness_ltable_cpp(phy)
    if (normalization == "tips" || normalization == TRUE) {
      n <- length(phy[, 1])
      betweenness_stat <- betweenness_stat / ((n - 2) * (n - 1))
    }
    return(betweenness_stat)
  }
  if (inherits(phy, "phylo")) {
    betweenness_stat <- calc_max_betweenness_cpp(phy)
    if (normalization == "tips" || normalization == TRUE) {
      n <- length(phy$tip.label)
      betweenness_stat <- betweenness_stat / ((n - 2) * (n - 1))
    }
    return(betweenness_stat)
  }
  stop("input object has to be phylo or ltable")
}
