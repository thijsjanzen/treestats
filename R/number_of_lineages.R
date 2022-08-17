#' Calculate the number of tips of a tree, including extinct tips.
#' @param phy phylo object
#' @return number of lineages
#' @export
number_of_lineages <- function(phy) {
  if (inherits(phy, "matrix")) {
    return(length(phy[, 1]))
  }

  if (inherits(phy, "phylo")) {
      return(length(phy$tip.label))
  }
  stop("input object has to be phylo or ltable")
}
