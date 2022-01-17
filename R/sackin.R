#' fast function using C++ to calculate the sackin index of (im)balance
#' @param phy phylogeny or ltable
#' @param normalization normalization, either 'none' (default), "yule" or "pda".
#' @return sackin index
#' @export
sackin <- function(phy, normalization = "none") {
  sackin_index <- apply_function_phy(phy, calc_sackin_cpp, normalization)
  return(sackin_index)
}
