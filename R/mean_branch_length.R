#' Calculates the mean branch length of a tree, including extinct branches.
#' @param phy phylo object or Ltable
#' @return mean branch length
#' @export
mean_branch_length <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_mean_branch_length_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_mean_branch_length_cpp(as.vector(phy$edge.length)))
  }

  stop("input object has to be phylo or ltable")
}
