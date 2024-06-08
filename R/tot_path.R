#' Total  path length
#' @description The total path length describes the sums of the depths
#' of all vertices of the tree.
#' @param phy phylo object or ltable
#' @return Total path length
#' @references  C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' @export
tot_path_length <- function(phy) {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    t_p <- sackin(phy) + tot_internal_path(phy)
    return(t_p)
  }
  stop("input object has to be phylo or ltable")
}
