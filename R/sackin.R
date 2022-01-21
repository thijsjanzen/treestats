#' fast function using C++ to calculate the sackin index of (im)balance
#' @param phy phylogeny or ltable
#' @param normalization normalization, either 'none' (default), "yule" or "pda".
#' @return sackin index
#' @references M. J. Sackin (1972). "Good" and "Bad" Phenograms.
#' Systematic Biology. 21:225-226.
#' @export
sackin <- function(phy, normalization = "none") {
  if (!ape::is.ultrametric(phy)) {
    stop("can only calculate sackin statistic for ultrametric tree")
  }
  sackin_index <- apply_function_phy(phy, calc_sackin_cpp, normalization)
  return(sackin_index)
}


