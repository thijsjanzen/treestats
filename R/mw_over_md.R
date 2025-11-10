#' Maximum width of branch depths divided by the maximum depth
#' @description Calculates the maximum width divided by the maximum depth.
#' @param phy phylogeny or ltable
#' @references C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
#' transmission patterns. Evolution, Medicine, and Public Health,
#' 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#' @export
mw_over_md <- function(phy) {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    return(calc_mw_over_md_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_mw_over_md_cpp(as.vector(t(phy$edge))))
  }

  stop("input object has to be phylo or ltable")
}
