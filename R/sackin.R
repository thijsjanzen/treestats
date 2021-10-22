#' fast function using C++ to calculate sackin index, please note that the
#' sackin index is not normalized here.
#' @param phy phylogeny or ltable
#' @return sackin index
#' @export
sackin <- function(phy) {
  sackin_index <- apply_function_ltab(phy, calc_sackin_cpp)
  return(sackin_index)
}
