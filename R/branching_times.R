#' Calculates branching times of a tree, using C++
#' @param phy phylo object or ltable
#' @return vector of branching times
#' @description C++ based alternative to `ape::branching.times`, please note
#' that to maximise speed, `treestats::branching_times` does not return
#' node names associated to the branching times, in contrast to the ape version.
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
