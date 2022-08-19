#' Fast function using C++ to calculate the diameter
#' @description The Diameter of a tree is defined as the maximum length of a
#' shortest path. When taking branch lengths into account, this is equal to
#' twice the crown age.
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @return Diameter
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." PloS one 16.12 (2021): e0259877.
#' @export
diameter <- function(phy,
                     weight = FALSE) {

  if (inherits(phy, "matrix")) {
    diam_stat <- calc_diameter_ltable_cpp(phy, weight)
    return(diam_stat)
  }
  if (inherits(phy, "phylo")) {
    diam_stat <- calc_diameter_cpp(phy, weight)
    return(diam_stat)
  }
  stop("input object has to be phylo or ltable")
}
