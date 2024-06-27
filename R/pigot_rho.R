#' Pigot's rho
#' @param phy phylo object
#' @return rho
#' @description Calculates the change in rate between the first half and the
#' second half of the extant phylogeny. Rho = (r2 - r1) / (r1 + r2), where r
#' reflects the rate in either the first or second half. The rate within a half
#' is given by (log(n2) - log(n1) / t, where n2 is the number of lineages at the
#' end of the half, and n1 the number of lineages at the start of the half. Rho
#' varies between -1 and 1, with a 0 indicating a constant rate across the
#' phylogeny, a rho < 0 indicating a slow down and a rho > 0 indicating a speed
#' up of speciation. In contrast to the Gamma statistic, Pigot's rho is not
#' sensitive to tree size.
#' @references Alex L. Pigot, Albert B. Phillimore, Ian P. F. Owens,
#' C. David L. Orme, The Shape and Temporal Dynamics of Phylogenetic Trees
#' Arising from Geographic Speciation, Systematic Biology, Volume 59, Issue 6,
#' December 2010, Pages 660–673, https://doi.org/10.1093/sysbio/syq058
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' pigot_rho(simulated_tree) # should be around 0.
#' ddd_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 7)$tes
#' pigot_rho(ddd_tree) # because of diversity dependence, should be < 0
pigot_rho <- function(phy) {

  if (inherits(phy, "matrix")) {
    return(calc_rho_ltable_cpp(phy))
  }

  if (inherits(phy, "phylo")) {
    if (!ape::is.ultrametric(phy)) {
      return(calc_rho_complete_cpp(phy))
    }
    return(calc_rho_cpp(phy))
  }

  stop("input object has to be phylo or ltable")
}
