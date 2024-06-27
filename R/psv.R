#' Phylogenetic Species Variability.
#' @description The phylogenetic species variability is bounded in [0, 1]. The
#' psv quantifies how phylogenetic relatedness decrease the variance of a
#' (neutral) trait shared by all species in the tree. As species become more
#' related, the psv tends to 0.
#' Please note that the psv is a special case of the Mean Pair Distance (see
#' appendix of Tucker et al. 2017 for a full derivation), and thus correlates
#' directly.
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", in which case the obtained mean
#' pairwise distance is normalized by the factor 2log(n), where n is the number
#' of tips.
#' @return Phylogenetic Species Variability
#' @references  Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007)
#' Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83
#'
#' Tucker, Caroline M., et al. "A guide to phylogenetic metrics for
#' conservation, community ecology and macroecology."
#' Biological Reviews 92.2 (2017): 698-715.
#' @export
psv <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {

    psv_stat <- calc_psv_cpp(as.vector(t(phy$edge)),
                                       phy$edge.length)

    if (normalization == "tips" || normalization == TRUE) {
      n <- length(phy$tip.label)
      psv_stat <- psv_stat / (2 * log(n))
    }
    return(psv_stat)
  }
  stop("input object has to be phylo or ltable")
}
