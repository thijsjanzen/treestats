#' Total path length
#' @description The total path length describes the sums of the depths
#' of all vertices of the tree. This can also be calculated as the sum of the
#' total internal path and the sackin index.
#' @param phy phylo object or ltable
#' @param normalization "none" or "yule". When providing "yule", the
#' statistic is normalized by the expectation under the Yule model, as given
#' in Bauman & Rosenberg 2026.
#' @return Total path length
#' @references  C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' Bauman, D. J., & Rosenberg, N. A. (2026). Basis statistics for phylogenetic
#' tree shape: tree balance formulas beyond expectations. arXiv preprint
#' @export
tot_path_length <- function(phy, normalization = "none") {

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)
  normalization <- check_normalization_key(normalization)
  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    t_p <- sackin(phy) + tot_internal_path(phy)

    if (normalization == "yule") {
      n <- length(phy$tip.label)
      norm_const <- 2 - 6 * n + 4 * n * harmonic(n)
      t_p <- t_p / norm_const
    }

    return(t_p)
  }
  stop("input object has to be phylo or ltable")
}
