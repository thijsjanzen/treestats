#' Maximum depth metric
#' @description The maximum depth metric, measures the maximal path (in edges),
#' between the tips and the root.
#' @param phy phylo object or ltable
#' @param normalization "none" or "tips", in which case the resulting statistic
#' is divided by the number of tips in the tree.
#' @return Maximum depth (in number of edges)
#' @references  C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' @export
max_depth <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    max_d_stat <- calc_max_depth_ltable_cpp(phy)
    if (normalization == "tips" || normalization == TRUE) {
      max_d_stat <- max_d_stat / length(phy[, 1])
    }
    return(max_d_stat)
  }
  if (inherits(phy, "phylo")) {
    max_d_stat <- calc_max_depth_cpp(as.vector(t(phy$edge)))
    if (normalization == "tips" || normalization == TRUE) {
      max_d_stat <- max_d_stat / length(phy$tip.label)
    }
    return(max_d_stat)
  }
  stop("input object has to be phylo or ltable")
}
