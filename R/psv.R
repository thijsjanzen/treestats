#' Fast function using C++ to calculate Phylogenetic Species Variability.
#' @description The phylogenetic species variability is bounded in [0, 1]. The
#' psv quantifies how phylogenetic relatedness decrease the variance of a
#' (neutral) trait shared by all species in the tree. As species become more
#' related, the psv tends to 0.
#' @param phy phylo object or ltable
#' @return Phylogenetic Species Variability
#' @references  Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007)
#' Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83
#' @export
psv <- function(phy) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {
    return(calc_psv_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}
