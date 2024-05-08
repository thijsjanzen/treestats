#' Height of a tree.
#' @description In a reconstructed tree, obtaining the tree height is fairly
#' straightforward, and the function beautier::get_crown_age does a great job
#' at it. However, in a non-ultrametric tree, that function no longer works.
#' Alternatively, taking the maximum value of adephylo::distRoot will also yield
#' the tree height (including the root branch), but will typically perform
#' many superfluous calculations and thus be slow.
#' @param phy phylo object
#' @return crown age
#' @export
tree_height <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(phy[1, 1])
  }

  if (inherits(phy, "phylo")) {
    crown_age_of_tree <- crown_age(phy)
    root_length <- 0
    if (!is.null(phy$root.edge)) {
      root_length <- phy$root.edge
    }
    tree_height <- crown_age_of_tree + root_length
    return(tree_height)
  }
  stop("input object has to be of class phylo")
}

#' Crown age of a tree.
#' @description In a reconstructed tree, obtaining the crown age is
#' fairly straightforward, and the function beautier::get_crown_age does
#' a great job at it. However, in a non-ultrametric tree, that function no
#' longer works. This function provides a functioning alternative
#' @param phy phylo object or ltable
#' @return crown age
#' @export
crown_age <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(phy[1, 1])
  }

  if (inherits(phy, "phylo")) {
    return(calc_crown_age_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
