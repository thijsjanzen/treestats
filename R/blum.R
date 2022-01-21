#' fast function using C++ to calculate the blum index of (im)balance
#' @param phy phylogeny or ltable
#' @return sackin index
#' @references M. G. B. Blum and O. Francois (2006). Which random processes
#' describe the Tree of Life? A large-scale study of phylogenetic tree
#' imbalance. Systematic Biology. 55:685-691.
#' @export
blum <- function(phy) {
  if (!ape::is.ultrametric(phy)) {
    stop("can only calculate sackin statistic for ultrametric tree")
  }
  blum_index <- apply_function_phy(phy, calc_blum_cpp)
  return(blum_index)
}


