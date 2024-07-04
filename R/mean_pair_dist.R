#' Mean Pairwise distance
#' @description Fast function using C++ to calculate the mean pairwise distance,
#' using the fast algorithm by Constantinos, Sandel & Cheliotis (2012).
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", in which case the obtained mean
#' pairwise distance is normalized by the factor 2log(n),
#' where n is the number of tips.
#' @return Mean pairwise distance
#' @references  Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and
#' Systematics 33:475-505.
#'
#' Tsirogiannis, Constantinos, Brody Sandel, and Dimitris Cheliotis.
#' "Efficient computation of popular phylogenetic tree measures." Algorithms in
#' Bioinformatics: 12th International Workshop, WABI 2012, Ljubljana, Slovenia,
#' September 10-12, 2012. Proceedings 12. Springer Berlin Heidelberg, 2012.
#' @export
mean_pair_dist <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }

  if (inherits(phy, "phylo")) {
    if (check_binary(phy)) {
      mpd <- calc_mpd_cpp(as.vector(t(phy$edge)),
                          phy$edge.length)
    } else {
      dist_mat <- ape::cophenetic.phylo(phy)
      diag(dist_mat) <- NA
      mpd <- mean(dist_mat, na.rm = TRUE)
    }
    if (normalization == "tips" || normalization == TRUE) {
      n <- length(phy$tip.label)
      mpd <- mpd / (2 * log(n))
    }
    return(mpd)
  }
  stop("input object has to be phylo or ltable")
}
