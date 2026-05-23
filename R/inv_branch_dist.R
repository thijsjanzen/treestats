#' Sum of Inverse branch length
#' @description This statistic calculates per tip, the sum of the inverse of
#' all branch lengths of the shortest path between the root and the tip.
#' @param phy phylo object or ltable
#' @param add_one should we take 1 / bl or 1 / (1 + bl) ? if TRUE, then
#' 1 / (1 + bl) is used. Default: FALSE.
#' @references Williams, P. H., Alonso-Alonso, P., Arbetman, M., Françoso,
#' E., Ghisbain, G., Huang, J., Orr, M. C., Ren, Z.-X., Streinzer, M., T
#' hanoosing, C., Vandame, R., Waite, M., & Brace, S. (2026). Evolutionary
#' Tree for All Bumblebee Species World-Wide Estimated by Combining
#' Information from Fast-Evolving Genes, Slow-Evolving Genes, and Genomic
#' Data (Apidae, Bombus). Insects, 17(6), 540.
#' https://doi.org/10.3390/insects17060540
#' @export
inv_branch_dist <- function(phy, add_one = FALSE) {

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = TRUE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    res <- calc_inv_path_cpp(as.vector(t(phy$edge)),
                             phy$edge.length,
                             add_one)
    distances <- res$distances
    names(distances) <- phy$tip.label[res$tip_ids]
    distances <- distances[sort(names(distances))]
    return(distances)
  }
  stop("input object has to be phylo or ltable")
}
