#' Area per pair index
#' @description The area per pair index calculates the sum of the number of
#' edges on the path between all two leaves. Instead, the area per pair index
#' (APP) can also be derived from the Sackin (S) and total cophenetic index
#' (TC):
#' \eqn{ APP = \frac{2}{n}\cdot S - \frac{4}{n(n-1)}\cdot TC}
#' \eqn{APP = 2/n * S - 4/(n(n-1)) * TC}
#' @param phy phylo object or ltable
#' @param normalization "none" or "yule", in which case the acquired result
#' is divided by the expectation for the Yule model.
#' @return Area per pair index
#' @references  T. Araújo Lima, F. M. D. Marquitti, and M. A. M. de Aguiar.
#' Measuring Tree Balance with Normalized Tree Area. arXiv e-prints, art.
#' arXiv:2008.12867, 2020.
#' @export
area_per_pair <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  n <- 0

  if (inherits(phy, "matrix")) {
    n <- length(phy[, 1])
  } else if (inherits(phy, "phylo")) {
    n <- length(phy$tip.label)
  } else {
    stop("input object has to be phylo or ltable")
  }

  result <- 2 / n * treestats::sackin(phy) -
            4 / (n * (n - 1)) * treestats::tot_coph(phy)

  if (normalization == "yule") {
    h_n <- sum(1 / (1:n))
    a <- (n + 1) / (n - 1)

    expected_value <- 4 * ((h_n - 1) * a - 1)
    result <- result / expected_value
  }

  return(result)
}
