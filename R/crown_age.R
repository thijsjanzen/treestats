#' Calculates the height of a tree using C++.
#' In a reconstructed tree, obtaining the tree height is fairly straightforward,
#' and the function beautier::get_crown_age does a great job at it.
#' However, in a non-ultrametric tree, that function no longer works.
#' Alternatively, taking the maximum value of adephylo::distRoot will also yield
#' the crown age, but will typically perform many superfluous calculations and
#' thus be slow.
#' @param phy phylo object
#' @return crown age
#' @export
tree_height <- function(phy) {
  crown_age_of_tree <- crown_age(phy)
  root_length <- 0
  if (!is.null(phy$root.edge)) {
    root_length <- phy$root.edge
  }
  tree_height <- crown_age_of_tree + root_length
  return(tree_height)
}

#' Calculates the height of a tree using C++.
#' In a reconstructed tree, obtaining the tree height is fairly straightforward,
#' and the function beautier::get_crown_age does a great job at it.
#' However, in a non-ultrametric tree, that function no longer works.
#' Alternatively, taking the maximum value of adephylo::distRoot will also yield
#' the crown age, but will typically perform many superfluous calculations and
#' thus be slow.
#' @param phy phylo object
#' @return crown age
#' @export
crown_age <- function(phy) {
  return(apply_function_phy(phy, calc_crown_age_cpp))
}

