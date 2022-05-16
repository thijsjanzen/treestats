#' Calculates the mean branch length of a tree, including extinct branches.
#' @param phy phylo object or Ltable
#' @return mean branch length
#' @export
mean_branch_length <- function(phy) {
  calc_mean_br_R <- function(focal_tree) {
    return(mean(focal_tree$edge.length, na.rm = TRUE))
  }

  return(calc_mean_br_R(phy))

  development <- "ongoing"

if (development == "done") {
  if (inherits(phy, "matrix")) {
    return(calc_mean_branch_length_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_mean_branch_length_cpp(phy))
  }
}
  stop("input object has to be phylo or ltable")
}
