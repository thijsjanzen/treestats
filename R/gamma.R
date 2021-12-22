#' @keywords internal
calc_gamma <- function(phy) {
  return(calc_gamma_cpp(phy))
}

#' calculate the gamma statistic, using C++.
#' @param phy phylo or multiPhylo object
#' @return gamma statistic
#' @export
gamma_statistic <- function(phy) {
  gamma_stat <- apply_function_phy(phy, calc_gamma)
  return(gamma_stat)
}
