#' Maximum width of branch depths.
#' @description Calculates the maximum width, this is calculated by first
#' collecting the depth of each node and tip across the entire tree, where the
#' depth represents the distance (in nodes) to the root. Then, the width
#' represents the number of occurrences of each possible depth. The maximal
#' width then returns the maximum number of such occurences.
#' @param phy phylogeny or ltable
#' @param normalization "none" or "tips", in which case the resulting statistic
#' is divided by the number of tips in the tree.
#' @return maximum width
#' @references C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' @export
max_width <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  if (inherits(phy, "matrix")) {
    max_w_stat <- calc_max_width_ltable_cpp(phy)
    if (normalization == "tips" || normalization == TRUE) {
      max_w_stat <- max_w_stat / length(phy[, 1])
    }
    return(max_w_stat)
  }
  if (inherits(phy, "phylo")) {
    max_w_stat <- calc_max_width_cpp(as.vector(t(phy$edge)))
    if (normalization == "tips" || normalization == TRUE) {
      max_w_stat <- max_w_stat / length(phy$tip.label)
    }
    return(max_w_stat)
  }

  stop("input object has to be phylo or ltable")
}
