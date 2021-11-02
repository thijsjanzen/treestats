#' @keywords internal
calc_pigot_rho_R <- function(phy) {

  ltab <- DDD::phylo2L(phy)
  crown_age <- ltab[1, 1]
  mid_point <- crown_age / 2

  n_crown <- 2
  n_extant <- sum(ltab[, 4] == -1)
  species_born_before_mid_point <- which(ltab[, 1] > mid_point)
  species_alive_after_mid_point <- which(ltab[, 4] == -1)
  species_mid_point <- species_born_before_mid_point %in% species_alive_after_mid_point
  num_species_mid_point <- sum(species_mid_point)

  r1 <- (log(num_species_mid_point) - log(n_crown)) / mid_point
  r2 <- (log(n_extant) - log(num_species_mid_point)) / mid_point

  rho <- (r2 - r1) / (r1 + r2)
  cat(n_crown, num_species_mid_point, n_extant, "\n")
  cat(r1, r2, rho, "\n")

  return(rho)
}

#' @keywords internal
calc_pigot_rho_cpp <- function(phy) {
  crown_age = max(ape::branching.times(phy))[[1]]
  return(calc_rho_cpp(phy, crown_age))
}

#' calculate Pigot's rho
#' @param phy phylo or multiPhylo object
#' @return rho
#' @description Calculates the change in rate between the first half and the
#' second half of the extant phylogeny. Rho = (r2 - r1) / (r1 + r2), where r
#' reflects the rate in either the first or second half. The rate within a half
#' is given by (log(n2) - log(n1) / t, where n2 is the number of lineages at the
#' end of the half, and n1 the number of lineages at the start of the half. Rho
#' varies between -1 and 1, with a 0 indicating a constant rate across the
#' phylogeny, a rho < 0 indicating a slow down and a rho > 0 indicating a speed
#' up of speciation. In contrast to the gamma statistic, Pigot's rho is not
#' sensitive to tree size.
#' @references Alex L. Pigot, Albert B. Phillimore, Ian P. F. Owens,
#' C. David L. Orme, The Shape and Temporal Dynamics of Phylogenetic Trees
#' Arising from Geographic Speciation, Systematic Biology, Volume 59, Issue 6,
#' December 2010, Pages 660–673, https://doi.org/10.1093/sysbio/syq058
#' @export
pigot_rho <- function(phy) {
  rho <- apply_function(phy, calc_pigot_rho_cpp)
  return(rho)
}
