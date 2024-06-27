#' Treeness statistic
#' @description Calculates the fraction of tree length on internal branches,
#' also known as treeness or stemminess
#' @param phy phylo object or Ltable
#' @return sum of all internal branch lengths (e.g. branches not leading
#' to a tip) divided by the sum over all branch lengths.
#' @export
treeness <- function(phy) {
  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {
    root_no <- min(phy$edge[, 1])

    int_branch_length <- sum(phy$edge.length[phy$edge[, 2] >= root_no])

    return(int_branch_length / (sum(phy$edge.length)))
  }

  stop("input object has to be phylo or ltable")
}
