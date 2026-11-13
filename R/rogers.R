#' Rogers J index of (im)balance.
#' @description The Rogers index is calculated as the total number of internal
#' nodes that are unbalanced, e.g. for which both daughter nodes lead to a
#' different number of extant tips. in other words, the number of nodes where
#' \eqn{L  !=  R} (where L(R) is the number of extant tips of the Left (Right)
#' daughter node).
#' @param phy phylo object or ltable
#' @param normalization "none", "yule", or "tips". When providing "yule", the
#' statistic is normalized by the expectation under the Yule model, as given
#' in Bauman & Rosenberg 2026. When providing "tips", the resulting statistic
#' is divided by the number of tips - 2 (e.g. the maximum value of the rogers
#' index for a tree).
#' @return Rogers index
#' @references  J. S. Rogers. Central Moments and Probability Distributions of
#' Three Measures of Phylogenetic Tree Imbalance. Systematic Biology,
#' 45(1):99-110, 1996. doi: 10.1093/sysbio/45.1.99.
#' Bauman, D. J., & Rosenberg, N. A. (2026). Basis statistics for phylogenetic
#' tree shape: tree balance formulas beyond expectations. arXiv preprint
#' arXiv:2610.40894.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#' unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#' rogers(balanced_tree)
#' rogers(unbalanced_tree) # should be higher
rogers <- function(phy, normalization = "none") {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  normalization <- check_normalization_key(normalization)

  normalize_yule <- function(stat, n) {
    norm_const <- 0
    if (n %% 2 == 0) { # even
        norm_const <- 2 * n * (harmonic((n / 2) - 1) - harmonic(n - 1) + 1)
    } else {  # uneven
        norm_const <- 2 * n * (harmonic((n - 1) / 2) - harmonic(n - 1) + 1) - 2
    }
    return(stat / norm_const)
  }


  if (inherits(phy, "matrix")) {
    rogers_stat <- calc_rogers_ltable_cpp(phy)
    if (normalization == "tips") {
      rogers_stat <- rogers_stat / (length(phy[, 1]) - 2)
    }
    if (normalization == "yule") {
      rogers_stat <- normalize_yule(rogers_stat, length(phy[, 1]))
    }
    return(rogers_stat)
  }
  if (inherits(phy, "phylo")) {
    rogers_stat <- calc_rogers_cpp(as.vector(t(phy$edge)))
    if (normalization == "tips") {
      rogers_stat <- rogers_stat / (length(phy$tip.label) - 2)
    }
    if (normalization == "yule") {
      rogers_stat <- normalize_yule(rogers_stat, length(phy$tip.label))
    }
    return(rogers_stat)
  }
  stop("input object has to be phylo or ltable")
}
