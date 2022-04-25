#' Fast function using C++ to calculate the mean pairwise distance.
#' @description The mean pairwise distance calculates the average distance
#' between all combinations of tips.
#' @param phy phylo object or ltable
#' @return Mean pairwise distance
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#' @export
mean_pair_dist <- function(phy) {
  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {
    return(calc_mpd_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
