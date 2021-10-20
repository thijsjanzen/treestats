
#' fast function using C++ to calculate sackin index
#' @param phy phylogeny or ltable
#' @return sackin index
#' @export
sackin <- function(phy) {
  sackin_index <- NA
  if (class(phy) == "multiPhylo") {
    phy <- lapply(phy, DDD::phylo2L)
    sackin_index <- lapply(phy, calc_sackin_cpp)
    sackin_index <- unlist(sackin_index)
  } else {

    if (!is_ltable(phy)) {
      phy <- DDD::phylo2L(phy)
    }

    sackin_index <- calc_sackin_cpp(phy)
  }
  return(sackin_index)
}
