#' Fast function using C++ to calculate the variance of leaf depth statistic
#' @description The variance of leaf depth statistic returns the variance
#' of depths across all tips.
#' @param phy phylo object or ltable
#' @return Variance of leaf depths
#' @references  T. M. Coronado, A. Mir, F. Rosselló, and L. Rotger.
#' On Sackin's original proposal: the variance of the leaves' depths as a
#' phylogenetic balance index. BMC Bioinformatics, 21(1), 2020.
#' doi: 10.1186/s12859-020-3405-1.
#' @export
var_leaf_depth <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_var_leaf_depth_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_var_leaf_depth_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}
