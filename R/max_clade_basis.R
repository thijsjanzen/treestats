#' Maximum Clade-basis statistic
#' @description The Maximum Clade-basis statistic measures consistency with the
#' Yule model, such that the statistic is maximized for a Yule tree. Values
#' lower than those expected under a Yule tree are indicative of the focal tree
#' not being the result of a Yule process. Although the statistic can also be
#' calculated for non-binary trees, the interpretation of the statistic is
#' unclear for non-binary trees.
#' @param phy phylogeny or ltable
#' @param normalization normalization, either 'none' (default) or yule
#' @return Statistic
#' @references Bauman, D. J., & Rosenberg, N. A. (2026). Basis statistics for
#' phylogenetic tree shape: tree balance formulas beyond expectations. arXiv
#' preprint arXiv:2610.40894.
#' @export
max_clade_basis <- function(phy, normalization = "none") {

  normalization <- check_normalization_key(normalization)

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    return(calc_max_clade_basis_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_max_clade_basis_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}
