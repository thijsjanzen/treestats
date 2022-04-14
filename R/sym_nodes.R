#' Fast function using C++ to calculate sym nodes metric
#' @description Balance metric that returns the total number of internal nodes
#' that are not-symmetric (confusingly enough). A node is considered symmetric
#' when both daughter trees have the same topology, measured as having the
#' same sum of depths, where depth is measured as the distance from the root
#' to the node/tip.
#' @param phy phylo object or ltable
#' @return Maximum depth (in number of edges)
#' @references  S. J. Kersting and M. Fischer. Measuring tree balance using
#' symmetry nodes — A new balance index and its extremal properties.
#' Mathematical Biosciences, page 108690, 2021. ISSN 0025-5564.
#' doi:https://doi.org/10.1016/j.mbs.2021.108690
#' @export
sym_nodes <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_var_leaf_depth_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {

    warning("calc_sym_nodes is unfinished, my sometimes give the wrong results")
    return(calc_sym_nodes_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}