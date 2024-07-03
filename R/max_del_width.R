#' Maximum difference of widths of a phylogenetic tree
#' @description Calculates the maximum difference of widths of a phylogenetic
#' tree. First, the widths are calculated by collecting the depth of each node
#' and tip across the entire tree, where the depth represents the distance
#' (in nodes) to the root. Then, the width represents the number of occurrences
#' of each possible depth. Then, we take the difference between each consecutive
#' width, starting with the first width. The maximum difference is then
#' returned - whereas the original statistic designed by Colijn and Gardy used
#' the absolute maximum difference, we here use the modified version as
#' introduced in Fischher 2023: this returns the maximum value, without
#' absoluting negative widths. This ensures that this metric is a proper
#' (im)balance metric, follwing Fischer 2023.
#' @param phy phylogeny or ltable
#' @param normalization "none" or "tips", in which case the resulting statistic
#' is divided by the number of tips in the tree.
#' @return maximum difference of widths
#' @references C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018..
#' Fischer, M., Herbst, L., Kersting, S., Kühn, A. L., & Wicke, K. (2023).
#' Tree Balance Indices: A Comprehensive Survey.
#' @export
max_del_width <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    max_dw_stat <- calc_max_del_width_ltable_cpp(phy)
    if (normalization == "tips" || normalization == TRUE) {
      max_dw_stat <- max_dw_stat / length(phy[, 1])
    }
    return(max_dw_stat)
  }
  if (inherits(phy, "phylo")) {
    max_dw_stat <- calc_max_del_width_cpp(as.vector(t(phy$edge)))
    if (normalization == "tips" || normalization == TRUE) {
      max_dw_stat <- max_dw_stat / length(phy$tip.label)
    }
    return(max_dw_stat)
  }

  stop("input object has to be phylo or ltable")
}
