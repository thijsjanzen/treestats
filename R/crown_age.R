#' Calculates the crown age of a tree using C++.
#' In a reconstructed tree, obtaining the crown age is fairly straightforward,
#' and the function beautier::get_crown_age does a great job at it.
#' However, in a non-ultrametric tree, that function no longer works.
#' Alternatively, taking the maximum value of adephylo::distRoot will also yield
#' the crown age, but will typically perform many superfluous calculations and
#' thus be slow.
#' @param phy phylo object
#' @return crown age
#' @export
crown_age <- function(phy) {

  if (is.null(phy$root.edge)) {
    phy$root.edge <- 0
  }

  return(apply_function_phy(phy, calc_crown_age_cpp))
}
