#' Max-bauman statistic of Yule-ness
#' @description The Maximum Bauman statistic measures consistency with the Yule
#' model, such that the statistic is maximized for a Yule tree. Values lower
#' than those expected under a Yule tree are indicative of the focal tree not
#' being the result of a Yule process
#' @param phy phylogeny or ltable
#' @param normalization normalization, either 'none' (default) or yule (not
#' implemented yet)
#' @return Max-bauman statistic
#' @references Bauman, D. J., & Rosenberg, N. A. (2026). Basis statistics for phylogenetic
#' tree shape: tree balance formulas beyond expectations. arXiv preprint
#' arXiv:2610.40894.
#' @export
max_bauman <- function(phy, normalization = "none") {

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    return(max_bauman_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}
