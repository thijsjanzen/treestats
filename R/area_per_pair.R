#' Fast function using C++ to calculate the area per pair index
#' @description The area per pair index calculates the sum of the number of
#' edges on the path between all two leaves. Instead, the area per pair index
#' can also be derived from the sackin and total cophenetic index.
#' @param phy phylo object or ltable
#' @return Total cophenetic index
#' @references  T. Araújo Lima, F. M. D. Marquitti, and M. A. M. de Aguiar.
#' Measuring Tree Balance with Normalized Tree Area. arXiv e-prints, art.
#' arXiv:2008.12867, 2020.
#' @export
area_per_pair <- function(phy) {
  n <- 0
  if (inherits(phy, "matrix")) {
    n <- length(phy[, 1])
  }
  if (inherits(phy, "phylo")) {
    n <- length(phy$tip.label)
  }
  result <- 2 / n * treestats::sackin(phy) -
            4 / (n * (n - 1)) * treestats::tot_coph(phy)
  return(result)
}
