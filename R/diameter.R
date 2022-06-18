#' Fast function using C++ to calculate the diameter
#' @description The Diameter of a tree is defined as the maximum length of a
#' shortest path. When taking branch lenghts into account, this is typically
#' twice the crown age. When not taking branch lengths into account, this may
#' be linked to max_depth
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @return Diameter
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
diameter <- function(phy, weight = TRUE) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    return(calc_diameter_cpp(phy, weight))
  }
  stop("input object has to be phylo or ltable")
}
