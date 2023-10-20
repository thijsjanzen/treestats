#' Fast function using C++ to calculate the total cophenetic index.
#' @description The total cophenetic index is the sum of the depth of the last
#' common ancestor of all pairs of leaves.
#' @param phy phylo object or ltable
#' @param normalization "none" or "yule", when "yule" is chosen, the statistic
#' is divided by the Yule expectation
#' @return Total cophenetic index
#' @references  A. Mir, F. Rosselló, and L. Rotger. A new balance index for
#' phylogenetic trees. Mathematical Bio-sciences, 241(1):125-136, 2013.
#' doi: 10.1016/j.mbs.2012.10.005.
#' @export
tot_coph <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    tot_coph_stat <- calc_tot_coph_ltable_cpp(phy)
    if (normalization == "yule" || normalization == TRUE) {
      n <- length(phy[, 1])
      if (n == 2) warning("normalization not valid for trees of size 2")
      h_n <- sum(1 / (1:n))
      yule_expected <-  n * (n + 1) - 2 * n * h_n
      tot_coph_stat <- tot_coph_stat / yule_expected
    }
    return(tot_coph_stat)
  }
  if (inherits(phy, "phylo")) {
    tot_coph_stat <- calc_tot_coph_cpp(as.vector(t(phy$edge)))
    if (normalization == "yule" || normalization == TRUE) {
      n <- length(phy$tip.label)
      if (n == 2) warning("normalization not valid for trees of size 2")
      h_n <- sum(1 / (1:n))
      yule_expected <-  n * (n + 1) - 2 * n * h_n
      tot_coph_stat <- tot_coph_stat / yule_expected
    }
    return(tot_coph_stat)
  }
  stop("input object has to be phylo or ltable")
}
