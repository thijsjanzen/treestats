#' Sum of Inverse branch length
#' @description This statistic calculates per tip, the sum of the inverse of
#' all branch lengths of the shortest path between the root and the tip.
#' @param phy phylo object or ltable
#' @param add_one should we take 1 / bl or 1 / (1 + bl) ? if TRUE, then
#' 1 / (1 + bl) is used. Default: TRUE
#' @return Distribution of the sum of inverse brnch length.
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
