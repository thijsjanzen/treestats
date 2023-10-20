#' Fast function using C++ to calculate the Wiener index
#' @description The Wiener index is defined as the sum of all shortest path
#' lengths between pairs of nodes in a tree.
#' @param phy phylo object or ltable
#' @param normalization if TRUE, the Wiener index is normalized by the number of
#' nodes, e.g. by choose(n, 2), where n is the number of nodes.
#' @param weight if TRUE, branch lenghts are used.
#' @return Wiener index
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' Mohar, B., Pisanski, T. How to compute the Wiener index of a graph.
#' J Math Chem 2, 267–277 (1988)
#' @export
wiener <- function(phy, normalization = FALSE, weight = TRUE) {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    return(calc_wiener_cpp(phy, normalization, weight))
  }
  stop("input object has to be phylo or ltable")
}
