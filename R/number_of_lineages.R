#' Calculate the number of tips of a tree, including extinct tips.
#' @param phy phylo object
#' @return number of lineages
#' @export
number_of_lineages <- function(phy) {

  calc_num_lin <- function(focal_tree) {
    return(focal_tree$Nnode + 1)
  }

  return(apply_function_phy(phy, calc_num_lin))
}
