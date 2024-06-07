#' Total internal path length
#' @description The total internal path length describes the sums of the depths
#' of all inner vertices of the tree.
#' @param phy phylo object or ltable
#' @return Total internal path length
#' @references  C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' @export
tot_internal_path <- function(phy) {

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE)

  if (inherits(phy, "matrix")) {
    max_d_stat <- calc_avg_depth_ltable_cpp(phy)
    return(max_d_stat)
  }
  if (inherits(phy, "phylo")) {
    t_i_p <- tot_internal_path_cpp(as.vector(t(phy$edge)))
    return(t_i_p)
  }
  stop("input object has to be phylo or ltable")
}
