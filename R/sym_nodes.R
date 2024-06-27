#' Symmetry nodes metric
#' @description Balance metric that returns the total number of internal nodes
#' that are not-symmetric (confusingly enough). A node is considered symmetric
#' when both daughter trees have the same topology, measured as having the
#' same sum of depths, where depth is measured as the distance from the root
#' to the node/tip.
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", in which case the resulting statistic
#' is divided by the number of tips - 2 (e.g. the maximum value of the symmetry
#' nodes index for a tree).
#' @return Maximum depth (in number of edges)
#' @references  S. J. Kersting and M. Fischer. Measuring tree balance using
#' symmetry nodes — A new balance index and its extremal properties.
#' Mathematical Biosciences, page 108690, 2021. ISSN 0025-5564.
#' doi:https://doi.org/10.1016/j.mbs.2021.108690
#' @export
sym_nodes <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    sym_nodes_stat <- calc_sym_nodes_cpp(as.vector(t(phy$edge)))
    if (normalization == "tips" || normalization == TRUE) {
      sym_nodes_stat <- sym_nodes_stat / (length(phy$tip.label) - 2)
    }
    return(sym_nodes_stat)
  }
  stop("input object has to be phylo or ltable")
}
