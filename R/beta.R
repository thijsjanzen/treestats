#' fast function using C++ to calculate the aldous beta statistic
#' @param phy phylogeny or ltable
#' @param upper_lim upper limit for beta parameter
#' @return sackin index
#' @export
beta_statistic <- function(phy, upper_lim = 10) {
  if (!ape::is.ultrametric(phy)) {
    stop("can only calculate beta statistic for ultrametric tree")
  }

  aldous_beta <- apply_function_phy(phy, calc_beta_cpp, upper_lim)
  return(aldous_beta)
}
