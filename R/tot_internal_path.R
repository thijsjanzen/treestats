#' Total internal path length
#' @description The total internal path length describes the sums of the depths
#' of all inner vertices of the tree.
#' @param phy phylo object or ltable
#' @return Total internal path length
#' @references  Knuth, Donald E. The Art of Computer Programming:
#' Fundamental Algorithms, volume 1. Addison-Wesley Professional, 1997.
#' @export
tot_internal_path <- function(phy) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    t_i_p <- tot_internal_path_cpp(as.vector(t(phy$edge)))
    return(t_i_p)
  }
  stop("input object has to be phylo or ltable")
}
