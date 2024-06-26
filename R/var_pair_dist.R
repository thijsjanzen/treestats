#' Fast function using C++ to calculate the variance of all pairwise distances.
#' @description After calculating all pairwise distances between all tips,
#' this function takes the variance across these values.
#' @param phy phylo object or ltable
#' @return Variance in pairwise distance
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#' @export
var_pair_dist <- function(phy) {

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE)

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
    return(calc_var_mpd_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
