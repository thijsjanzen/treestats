#' Intensive quadratic entropy statistic J.
#' @description The intensive quadratic entropy statistic J is given by the
#' average distance between two randomly chosen species, thus given by the
#' sum of all pairwise distances, divided by \eqn{S^2}, where S is the number
#' of tips of the tree.
#' @param phy phylo object or ltable
#' @return intensive quadratic entropy statistic J
#' @references  Izsák, János, and Laszlo Papp. "A link between ecological
#' diversity indices and measures of biodiversity."
#' Ecological Modelling 130.1-3 (2000): 151-156.
#' @export
entropy_j <- function(phy) {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {
    if (ape::is.rooted(phy)) {
       return(calc_J_cpp(as.vector(t(phy$edge)),
                      phy$edge.length))
    } else {
      n <- length(phy$tip.label)
      return(treestats::mean_pair_dist(phy) / n)
      # notice that the statistic divides by S^2, but since the mean is
      # already dividing by S, this is fine.
    }
  }
  stop("input object has to be phylo or ltable")
}
