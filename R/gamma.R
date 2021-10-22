#' gamma statistic, from ape package
#' @param phy phylo or multiPhylo object
#' @return gamma statistic
#' @export
gamma_statistic <- function(phy) {
  gamma_stat <- apply_function(phy, ape::gammaStat)
  return(gamma_stat)
}
