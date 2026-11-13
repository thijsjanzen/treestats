#' Total internal path length
#' @description The total internal path length describes the sums of the depths
#' of all inner vertices of the tree.
#' @param normalization "none" or "yule". When providing "yule", the
#' statistic is normalized by the expectation under the Yule model, as given
#' in Bauman & Rosenberg 2026.
#' @param phy phylo object or ltable
#' @return Total internal path length
#' @references  Knuth, Donald E. The Art of Computer Programming:
#' Fundamental Algorithms, volume 1. Addison-Wesley Professional, 1997.
#' Bauman, D. J., & Rosenberg, N. A. (2026). Basis statistics for phylogenetic
#' tree shape: tree balance formulas beyond expectations. arXiv preprint
#' arXiv:2610.40894.
#' @export
tot_internal_path <- function(phy, normalization = "none") {

  normalization <- check_normalization_key(normalization)
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    t_i_p <- tot_internal_path_cpp(as.vector(t(phy$edge)))

    if (normalization == "yule") {
      n <- length(phy$tip.label)
      norm_const <- 2 - 4 * n + 2 * n * harmonic(n)
      t_i_p <- t_i_p / norm_const
    }

    return(t_i_p)
  }
  stop("input object has to be phylo or ltable")
}
