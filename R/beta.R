#' fast function using C++ to calculate the aldous beta statistic
#' @param phy phylogeny or ltable
#' @param upper_lim upper limit for beta parameter
#' @param algorithm optimization algorithm used, default is "subplex", also
#' available are "COBYLA" (Constrained Optimization BY Linear Approximations)
#' and "simplex".
#' @param abs_tol absolute stopping criterion. Default is 1e-4.
#' @param rel_tol relative stopping criterion. Default is 1e-6
#' @return sackin index
#' @export
beta_statistic <- function(phy,
                           upper_lim = 10,
                           algorithm = "subplex",
                           abs_tol = 1e-4,
                           rel_tol = 1e-6) {
  if (!ape::is.ultrametric(phy)) {
    stop("can only calculate beta statistic for ultrametric tree")
  }

  aldous_beta <- apply_function_phy(phy, calc_beta_cpp, upper_lim,
                                    algorithm, abs_tol, rel_tol)
  return(aldous_beta)
}
