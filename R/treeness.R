#' Treeness statistic
#' @description Calculates the fraction of tree length on internal branches,
#' also known as treeness or stemminess
#' @param phy phylo object or Ltable
#' @return sum of all internal branch lengths (e.g. branches not leading
#' to a tip) divided by the sum over all branch lengths.
#' @export
treeness <- function(phy) {

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    # in an unrooted tree, this number still indicates the first non-tip
    int_branch_length <- get_internal_branches(phy)
    return(sum(int_branch_length) / (sum(phy$edge.length)))
  }

  stop("input object has to be phylo or ltable")
}
