#' Calculates branching times of a tree, using C++
#' @param phy phylo object or ltable
#' @return vector of branching times
#' @export
branching_times <- function(phy) {
  apply_function_phy_ltable(phy,
                            branching_times_cpp,
                            branching_times_ltable_cpp,
                            only_extant = FALSE)
}
