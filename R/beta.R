
#' fast function using C++ to calculate the aldous beta statistic
#' @param phy phylogeny or ltable
#' @param upper_lim upper limit for beta parameter
#' @return sackin index
#' @export
beta_statistic <- function(phy, upper_lim = 10) {
  aldous_beta <- apply_function_ltab(phy, calc_beta_cpp, upper_lim)
  return(aldous_beta)
}
