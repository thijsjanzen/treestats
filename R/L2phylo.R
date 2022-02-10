#' Convert an L table to phylo object
#' @param ltab ltable
#' @param drop_extinct should extinct species be dropped from the phylogeny?
#' @return phylo object
#' @export
l_to_phylo <- function(ltab, drop_extinct = TRUE) {

  newick_str <- l_to_newick(ltab, drop_extinct)
  phylo_tree <- ape::read.tree(text = newick_str)

  return(phylo_tree)
}

#' Convert an L table to newick string
#' @param ltab ltable
#' @param drop_extinct should extinct species be dropped from the phylogeny?
#' @return phylo object
#' @export
ltable_to_newick <- function(ltab, drop_extinct = TRUE) {
  newick_str <- l_to_newick(ltab, drop_extinct)
  return(newick_str)
}
