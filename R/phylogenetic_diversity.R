#' @keywords internal
calc_phylogenetic_diversity <- function(phy, t, extinct_tol) {
  if (t == 0 && ape::is.ultrametric(phy, option = 2)) {
    return(sum(phy$edge.length)) # no need to pass to Rcpp
  } else {
    if (ape::is.rooted(phy)) {
      return(calc_phylodiv_cpp(phy, t, extinct_tol))
    }
  }
}



#' Phylogenetic diversity at time point t
#' @description The phylogenetic diversity at time t is given by the total
#' branch length of the tree reconstructed up until time point t. Time is
#' measured increasingly, with the crown age equal to 0. Thus, the time at
#' the present is equal to the crown age.
#' @param input_obj phylo object or Ltable
#' @param t time point at which to measure phylogenetic diversity, alternatively
#' a vector of time points can also be provided. Time is measured with 0 being
#' the present.
#' @param extinct_tol tolerance to determine if a lineage is extinct at time t.
#' Default is 1/100 * smallest branch length of the tree.
#' @return phylogenetic diversity, or vector of phylogenetic diversity measures
#' if a vector of time points is used as input.
#' @references Faith, Daniel P. "Conservation evaluation and phylogenetic
#' diversity." Biological conservation 61.1 (1992): 1-10.
#' @export
phylogenetic_diversity <- function(input_obj,
                                   t = 0,
                                   extinct_tol = NULL) {

  if (inherits(input_obj, "matrix")) {

    if (length(t) == 1) {
      if (t != 0) {
        stop("Ltable implemenation can only be used for t = present = 0")
      }
      if (sum(input_obj[, 4] > 0)) {
        stop("Ltable implementation only supports extant trees")
      }

      return(calc_phylodiv_ltable_cpp(input_obj))
    } else {
       stop("Ltable implemenation can only be used for a single time point, t = 0") #nolint
    }
  }

  if (inherits(input_obj, "phylo")) {
    if (!ape::is.rooted(input_obj)) {
      stop("phylogenetic diversity is not available for unrooted trees")
    }

    if (is.null(extinct_tol)) {
      extinct_tol <- min(input_obj$edge.length) / 100
    }

    if (length(t) == 1) {
      return(calc_phylogenetic_diversity(input_obj, t, extinct_tol))
    }

    fun_to_apply <- function(focal_time) {
      return(calc_phylogenetic_diversity(input_obj,
                                         focal_time,
                                         extinct_tol))
    }

    out <- lapply(t, fun_to_apply)
    return(unlist(out))
  }

  stop("input object has to be phylo or ltable")
}
