#' gamma statistic, from ape package
#' @param phy phylo or multiPhylo object
#' @return gamma statistic
#' @export
gamma <- function(phy) {
  gamma <- apply_function(phy, ape::gammaStat)
  return(gamma)
}
