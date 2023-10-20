#' Fast function using C++ to calculate the B1 metric
#' @description Balance metric (in the case of a binary tree), which measures
#' the sum across all internal nodes of one over the maximum depth of all
#' attached tips to that node.
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", in which case the resulting
#' statistic is divided by the number of tips in the tree, as a crude way of
#' normalization.
#' @return B1 statistic
#' @references  K.-T. Shao and R. R. Sokal. Tree Balance.
#' Systematic Zoology, 39(3):266, 1990. doi: 10.2307/2992186.
#' @export
b1 <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)
  if (inherits(phy, "matrix")) {
    b1_stat <- calc_b1_ltable_cpp(phy)
    if (normalization == "tips" || normalization == TRUE) {
      n <- length(phy[, 1])
      b1_stat <- b1_stat / n
    }
    return(b1_stat)
  }
  if (inherits(phy, "phylo")) {
    b1_stat <- calc_b1_cpp(as.vector(t(phy$edge)))
    if (normalization == "tips" || normalization == TRUE) {
      n <- length(phy$tip.label)
      b1_stat <- b1_stat / n
    }
    return(b1_stat)
  }
  stop("input object has to be phylo or ltable")
}
