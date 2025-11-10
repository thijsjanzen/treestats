#' Average vertex depth metric
#' @description The average vertex depth metric, measures the average path
#' (in edges), between the tips and the root.
#' @param phy phylo object or ltable
#' @return Average depth (in number of edges)
#' @references  C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' @export
avg_vert_depth <- function(phy) {
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)


  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    max_d_stat <- calc_avg_vert_depth_cpp(as.vector(t(phy$edge)))
    return(max_d_stat)
  }
  stop("input object has to be phylo or ltable")
}
