#' Calculates the mean branch length of a tree, including extinct branches.
#' @param phy phylo object
#' @return mean branch length
#' @export
mean_branch_length <- function(phy) {

  calc_mean_br <- function(focal_tree) {
    return(mean(focal_tree$edge.length, na.rm = TRUE))
  }

  return(apply_function_phy(phy, calc_mean_br))
}
