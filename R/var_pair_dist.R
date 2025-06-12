#' Variance of all pairwise distances.
#' @description After calculating all pairwise distances between all tips,
#' this function takes the variance across these values.
#' @param phy phylo object or ltable
#' @return Variance in pairwise distance
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#' @export
var_pair_dist <- function(phy) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy)
  }
  if (inherits(phy, "phylo")) {

    if (ape::is.rooted(phy)) {
       return(calc_var_mpd_cpp(phy))
    } else {
      dist_mat <- ape::cophenetic.phylo(phy)
      dist_mat <- dist_mat[lower.tri(dist_mat)]
      n <- length(dist_mat)
      var_mpd <- stats::var(dist_mat, na.rm = TRUE, use = "everything")
      # var uses sample variance, we use population variance
      var_mpd <- var_mpd * (n - 1) / n
      return(var_mpd)
    }
  }
  stop("input object has to be phylo or ltable")
}
