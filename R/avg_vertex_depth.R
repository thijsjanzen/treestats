#' Average vertex depth metric
#' @description The average vertex depth metric, measures the average path
#' (in edges), between the tips and the root.
#' @param phy phylo object or ltable
#' @param normalization "none" or "yule". When providing "yule", the
#' statistic is normalized by the expectation under the Yule model, as given
#' in Bauman & Rosenberg 2026.
#' @return Average depth (in number of edges)
#' @references  C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' Bauman, D. J., & Rosenberg, N. A. (2026). Basis statistics for phylogenetic
#' tree shape: tree balance formulas beyond expectations. arXiv preprint
#' arXiv:2610.40894.
#' @export
avg_vert_depth <- function(phy, normalization = "none") {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    # polytomies are checked at the Rcpp side
    max_d_stat <-  calc_avg_vert_depth_cpp(as.vector(t(phy$edge)))
    if (normalization == "yule") {
      n <- length(phy$tip.label)
      norm_const <- (2 - 6 * n + 4 * n * harmonic(n)) / (2 * n - 1)
      max_d_stat <- max_d_stat / norm_const
    }
    return(max_d_stat)
  }
  stop("input object has to be phylo or ltable")
}
