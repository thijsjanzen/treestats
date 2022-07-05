#' Fast function using C++ to calculate the diameter
#' @description The Diameter of a tree is defined as the maximum length of a
#' shortest path. When taking branch lengths into account, this is equal to
#' twice the crown age.
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @param normalization "none" or "minmax", where the found value is rescaled
#' using the minimum and maximum diameter value expected for a tree of n tips,
#' where the minimum is n and the maximum is given by 2log_2(n).
#' @return Diameter
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
diameter <- function(phy, weight = FALSE, normalization = "none") {

  if (inherits(phy, "matrix")) {
    diam_stat <- calc_diameter_ltable_cpp(phy, weight)
    if (normalization == "minmax") {
      n <- length(phy[, 1])
      min_val <- n
      max_val <- 2 * log2(n)
      diam_stat <- (diam_stat - min_val) / (max_val - min_val)
    }
    return(diam_stat)
  }
  if (inherits(phy, "phylo")) {
    diam_stat <- calc_diameter_cpp(phy, weight)
    if (normalization == "minmax") {
      n <- length(phy$tip.label)
      min_val <- n
      max_val <- 2 * log2(n)
      diam_stat <- (diam_stat - min_val) / (max_val - min_val)
    }
    return(diam_stat)
  }
  stop("input object has to be phylo or ltable")
}
