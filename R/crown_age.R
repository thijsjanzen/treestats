#' Calculates the crown age of a tree.
#' @param phy phylo object
#' @return crown age
#' @export
crown_age <- function(phy) {
  return(apply_function_phy(phy, calc_crown_age_cpp))
}
