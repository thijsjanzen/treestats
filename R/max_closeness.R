#' Fast function using C++ to calculate maximum closeness
#' @description Closeness is defined as 1 / Farness, where Farness is the sum
#' of distances from a node to all the other nodes in the tree. Here, we return
#' the node with maximum closeness.
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @return Maximum Closeness
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' Wang W, Tang CY. Distributed computation of classic and exponential closeness
#' on tree graphs. Proceedings of the American Control Conference. IEEE; 2014.
#' p. 2090–2095.
#' @export
max_closeness <- function(phy, weight = TRUE) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    return(calc_max_closeness_cpp(phy, weight))
  }
  stop("input object has to be phylo or ltable")
}
