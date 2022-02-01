#' Calculates phylogenetic diversity at time point t, e.g. the total branch
#' length of the tree reconstructed up until time point t. Time is measured
#' increasingly, with the crown age equal to 0. Thus, the time at the present
#' is equal to the crown age.
#' @param phy phylo object
#' @param t time point at which to measure phylogenetic diversity, alternatively
#' a vector of time points can also be provided.
#' @param extinct_tol tolerance to determine if a lineage is extinct at time t.
#' Default is 1/100 * smallest branch length of the tree.
#' @return phylogenetic diversity, or vector of phylogenetic diversity measures
#' if a vector of time points is used as input.
#' @export
phylogenetic_diversity <- function(phy,
                                   t = 1e10,
                                   extinct_tol = NULL) {
  crown_age <- max(treestats::branching_times(phy))[[1]]
  if (is.null(extinct_tol)) {
    extinct_tol <- min(phy$edge.length) / 100
  }

  if (length(t) == 1) {
    return(calc_phylodiv_cpp(phy, t, crown_age, extinct_tol))
  }


  fun_to_apply <- function(focal_time) {
    return(calc_phylodiv_cpp(phy, focal_time, crown_age, extinct_tol))
  }

  out <- lapply(t, fun_to_apply)
  return(unlist(out))
}
