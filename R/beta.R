
#' fast function using C++ to calculate the aldous beta statistic
#' @param phy phylogeny or ltable
#' @return sackin index
#' @export
beta <- function(phy) {
  aldous_beta <- NA
  if (class(phy) == "multiPhylo") {
    phy <- lapply(phy, DDD::phylo2L)
    aldous_beta <- lapply(phy, calc_beta_cpp)
    aldous_beta <- unlist(aldous_beta)
  } else {

    if (!is_ltable(phy)) {
      phy <- DDD::phylo2L(phy)
    }

    aldous_beta <- calc_beta_cpp(phy)
  }
  return(aldous_beta)
}
