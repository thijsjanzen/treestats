#' calculate the mean branch length
#' @param phy phylo or multiPhylo object
#' @return mean branch length
#' @export
mean_branch_length <- function(phy) {

  calc_mean_br <- function(focal_tree) {
    return(mean(focal_tree$edge.length, na.rm = TRUE))
  }

  return(apply_function(phy, calc_mean_br))
}
