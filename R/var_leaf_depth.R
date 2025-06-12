#' Variance of leaf depth statistic
#' @description The variance of leaf depth statistic returns the variance
#' of depths across all tips, where depth of a tip indicates the distance of
#' the tip to the root.
#' @param phy phylo object or ltable
#' @param normalization "none" or "yule", when "yule" is chosen, the statistic
#' is divided by the Yule expectation
#' @return Variance of leaf depths
#' @references  T. M. Coronado, A. Mir, F. Rosselló, and L. Rotger.
#' On Sackin's original proposal: the variance of the leaves' depths as a
#' phylogenetic balance index. BMC Bioinformatics, 21(1), 2020.
#' doi: 10.1186/s12859-020-3405-1.
#' @export
var_leaf_depth <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    var_leaf_depth_stat <- calc_var_leaf_depth_ltable_cpp(phy)
    if (normalization == "yule" || normalization == TRUE) {
      n <- length(phy[, 1])
      yule_expected <-  2 * log(n)
      var_leaf_depth_stat <- var_leaf_depth_stat / yule_expected
    }
    return(var_leaf_depth_stat)
  }
  if (inherits(phy, "phylo")) {
    var_leaf_depth_stat <- calc_var_leaf_depth_cpp(as.vector(t(phy$edge)))
    if (normalization == "yule" || normalization == TRUE) {
      n <- length(phy$tip.label)
      yule_expected <-  2 * log(n)
      var_leaf_depth_stat <- var_leaf_depth_stat / yule_expected
    }
    return(var_leaf_depth_stat)
  }
  stop("input object has to be phylo or ltable")
}
