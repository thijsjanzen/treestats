#' Mean I statistic.
#' @description The mean I value is defined for all nodes with at least 4 tips
#' connected, such that different topologies can be formed. Then, for each node,
#' \eqn{I = (nm - nt/2) / (nt - 1 - nt/2)}, where nt is the total number of tips
#' descending from that node, nm is the daughter branch leading to most tips,
#' and \eqn{nt/2} is the minimum size of the maximum branch, rounded up.
#' Following Purvis et al 2002, we perform a correction on I, where we correct
#' I for odd nt, such that \eqn{I' = I * (nt - 1) / nt}. This correction ensures
#' that I is independent of nt. We report the mean value across all I' (again,
#' following Purvis et al. 2002).
#' @param phy phylo object or ltable
#' @return average I value across all nodes
#' @references  G. Fusco and Q. C. Cronk. A new method for evaluating the shape
#' of large phylogenies. Journal of Theoretical Biology, 1995.
#' doi: 10.1006/jtbi.1995.0136.
#'
#' A. Purvis, A. Katzourakis, and P.-M. Agapow. Evaluating Phylogenetic Tree
#' Shape: Two Modifications to Fusco & Cronks Method. Journal of Theoretical
#' Biology, 2002. doi: 10.1006/jtbi.2001.2443.
#' @export
mean_i <- function(phy) {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    if (length(phy[, 1]) < 4) {
      warning("I statistic is only available for trees with at least 4 tips.")
      return(NA)
    }
    return(calc_Ibased_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {

    if (length(phy$tip.label) < 4) {
      warning("I statistic is only available for trees with at least 4 tips.")
      return(NA)
    }

    return(calc_Ibased_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}


#' Mean I statistic.
#' @description The mean I value is defined for all nodes with at least 4 tips
#' connected, such that different topologies can be formed. Then, for each node,
#' \eqn{I = (nm - nt/2) / (nt - 1 - nt/2)}, where nt is the total number of tips
#' descending from that node, nm is the daughter branch leading to most tips,
#' and \eqn{nt/2} is the minimum size of the maximum branch, rounded up.
#' Following Purvis et al 2002, we perform a correction on I, where we correct
#' I for odd nt, such that \eqn{I' = I * (nt - 1) / nt}. This correction ensures
#' that I is independent of nt. We report the mean value across all I' (again,
#' following Purvis et al. 2002).
#' @param phy phylo object or ltable
#' @return average I value across all nodes
#' @references  G. Fusco and Q. C. Cronk. A new method for evaluating the shape
#' of large phylogenies. Journal of Theoretical Biology, 1995.
#' doi: 10.1006/jtbi.1995.0136.
#'
#' A. Purvis, A. Katzourakis, and P.-M. Agapow. Evaluating Phylogenetic Tree
#' Shape: Two Modifications to Fusco & Cronks Method. Journal of Theoretical
#' Biology, 2002. doi: 10.1006/jtbi.2001.2443.
#' @export
i_stat <- function(phy) {
  return(mean_i(phy))
}
