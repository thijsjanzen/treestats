#' Fast function using C++ to calculate maximum betweenness centrality.
#' @description Betweenness centrality associates with each node v, the two
#' nodes u, w, for which the shortest path between u and w runs through v, if
#' the tree were re-rooted at node v. Then, we report the node with maximum
#' betweenness centrality.
#' @param phy phylo object or ltable
#' @return Maximum Betweenness
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
max_betweenness <- function(phy) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {
    return(calc_max_betweenness_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
