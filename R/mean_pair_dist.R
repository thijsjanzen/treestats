#' Fast function using C++ to calculate the mean pairwise distance.
#' @description The mean pairwise distance calculates the average distance
#' between all combinations of tips.
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", in which case the obtained mean
#' pairwise distance is normalized by the factor 2log(n), where n is the number
#' of tips.
#' @return Mean pairwise distance
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#' @export
mean_pair_dist <- function(phy, normalization = "none") {
  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {
    n <- length(phy$tip.label)
    m <- phy$Nnode
    nm <- n + m
    if (nm > 46340) { # sqrt(2^31 - 1) #nolint
      stop("tree too big")
    }
    mpd <- calc_mpd_cpp(phy)
    if (normalization == "tips") {
      n <- length(phy$tip.label)
      mpd <- mpd / (2 * log(n))
    }
    return(mpd)
  }
  stop("input object has to be phylo or ltable")
}
