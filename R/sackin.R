#' fast function using C++ to calculate sackin index, please note that the
#' sackin index is not normalized here.
#' @param phy phylogeny or ltable
#' @param norm normalization, either 'none' (default), "yule" or "pda".
#' @return sackin index
#' @export
sackin <- function(phy, norm = 'none') {
  sackin_index <- apply_function_ltab(phy, calc_sackin_cpp, norm)
  return(sackin_index)
}
