#' Diameter statistic
#' @description The Diameter of a tree is defined as the maximum length of a
#' shortest path. When taking branch lengths into account, this is equal to
#' twice the crown age.
#' When the tree is unrooted, we add 1 to the unweighted diameter, to reflect
#' traversing the (virtual) root.
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @return Diameter
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." PloS one 16.12 (2021): e0259877.
#' @export
diameter <- function(phy,
                     weight = FALSE) {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  if (inherits(phy, "matrix")) {
    diam_stat <- calc_diameter_ltable_cpp(phy, weight)
    return(diam_stat)
  }
  if (inherits(phy, "phylo")) {
    diam_stat <- NA
    if (ape::is.rooted(phy)) {
      diam_stat <- calc_diameter_cpp(phy, weight)
    } else {
      if (!weight) {
        phy$edge.length <- rep(1, length(phy$edge.length))
      }
      dist_mat <- ape::cophenetic.phylo(phy)
      diam_stat <- max(dist_mat, na.rm = TRUE)
    }
    return(diam_stat)
  }
  stop("input object has to be phylo or ltable")
}
