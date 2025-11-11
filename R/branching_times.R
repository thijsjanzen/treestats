#' Branching times of a tree
#' @param phy phylo object or ltable
#' @return vector of branching times
#' @description C++ based alternative to
#' \link[ape:branching.times]{ape::branching.times}.
#' @export
branching_times <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(branching_times_ltable_cpp(phy))
  }

  if (inherits(phy, "phylo")) {
    answ <- branching_times_cpp(phy)
    n <- length(phy$tip.label)
    if (is.null(phy$node.label)) {
       names(answ) <- (n + 1):(n + phy$Nnode)
    } else {
       names(answ) <- phy$node.label
    }
    return(answ)
  }
  stop("input object has to be phylo or ltable")
}
