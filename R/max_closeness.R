#' Maximum closeness
#' @description Closeness is defined as 1 / Farness, where Farness is the sum
#' of distances from a node to all the other nodes in the tree. Here, we return
#' the node with maximum closeness.
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @param normalization "none" or "tips", in which case an arbitrary post-hoc
#' correction is performed by dividing by the expectation of n log(n), where n
#' is the number of tips.
#' @return Maximum Closeness
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' Wang W, Tang CY. Distributed computation of classic and exponential closeness
#' on tree graphs. Proceedings of the American Control Conference. IEEE; 2014.
#' p. 2090–2095.
#' @export
max_closeness <- function(phy, weight = TRUE, normalization = "none") {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {

    if (ape::is.binary(phy) && ape::is.rooted(phy)) {
      closeness_stat <- calc_max_closeness_cpp(phy, weight)
    } else {
      if (!weight) phy$edge.length <- rep(1, length(phy$edge.length))
      node_dist <- ape::dist.nodes(phy)
      closeness_stat <- max(1 / rowSums(node_dist))
    }

    if (normalization == "tips" || normalization == TRUE) {
      n <- length(phy$edge)
      expectation <- 1.0 / (n * log(n))
      closeness_stat <- closeness_stat / expectation
    }
    return(closeness_stat)
  }
  stop("input object has to be phylo or ltable")
}
