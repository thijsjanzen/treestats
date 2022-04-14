#' Fast function using C++ to calculate the B2 metric
#' @description Balance metric that uses the Shannon-Wiener statistic of
#' information content. The b2 measure is given by the sum over the depths of
#' all tips, divided by 2^depth: sum Ni / 2^Ni
#' @param phy phylo object or ltable
#' @return Maximum depth (in number of edges)
#' @references  K.-T. Shao and R. R. Sokal. Tree Balance.
#' Systematic Zoology, 39(3):266, 1990. doi: 10.2307/2992186.
#' @export
b2 <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_b2_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_b2_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}
