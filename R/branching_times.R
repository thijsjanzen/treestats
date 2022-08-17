#' Calculates branching times of a tree, using C++
#' @param phy phylo object or ltable
#' @return vector of branching times
#' @export
branching_times <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(branching_times_ltable_cpp(phy))
  }

  if (inherits(phy, "phylo")) {
    return(branching_times_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
