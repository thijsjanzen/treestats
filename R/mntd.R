#' Mean Nearest Taxon distance
#' @description Per tip, evaluates the shortest distance to another tip,
#' then takes the average across all tips.
#' @param phy phylo object or ltable
#' @return Mean Nearest Taxon Distance.
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#' @export
mntd <- function(phy) {
  # only defined for extant trees!
  if (inherits(phy, "matrix")) {
    if (sum(phy[, 4] != -1)) {
      # we now have to do it brute force:
      phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
      return(brute_mntd(phy))
    }
    return(calc_mntd_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    if (!ape::is.ultrametric(phy, option = 2)) {
      return(brute_mntd(phy)) # slower though
    }
    return(calc_mntd_cpp(phy))
  }
  stop("input object has to be phylo or ltable")
}

#' @keywords internal
brute_mntd <- function(phy) {
  dist_mat <- ape::cophenetic.phylo(phy)
  diag(dist_mat) <- NA
  ntd <- apply(dist_mat, 1,
               function(x) {
                 return(min(x, na.rm = TRUE))
               }
              )
  return(mean(ntd))
}
