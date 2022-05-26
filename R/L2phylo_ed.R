#' Convert an reverse-timescale L table to phylo object,
#' only use it with PDD simulation
#' @param ltab ltable
#' @param t simulation time when converting it to phylogeny
#' @param drop_extinct should extinct species be dropped from the phylogeny?
#' @return phylo object
#' @export
l_to_phylo_ed <- function(ltab, t, drop_extinct = TRUE) {

  newick_str <- l_to_newick_ed(ltab, t, drop_extinct)
  phylo_tree <- ape::read.tree(text = newick_str)

  return(phylo_tree)
}

#' Convert an reverse-timescale L table to newick format string
#' @param ltab ltable
#' @param t simulation time when converting it to phylogeny
#' @param drop_extinct should extinct species be dropped from the phylogeny?
#' @return newick string
#' @export
l_to_newick_ed <- function(ltab, t, drop_extinct = TRUE) {

  newick_str <- l_to_newick_ed(ltab, t, drop_extinct)

  return(newick_str)
}
