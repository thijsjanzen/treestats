#' Blum index of (im)balance.
#' @description The Blum index of imbalance (also known as the s-shape
#' statistic (see \link{sshape})) calculates the sum of \eqn{log(N-1)} over all
#' internal nodes, where N represents the total number of extant tips connected
#' to that node. An alternative implementation can be found in the Castor R
#' package.
#' @param phy phylogeny or ltable
#' @param normalization "none" or "yule", in which case the acquired result
#' is divided by the expectation for the Yule model.
#' @return Blum index of imbalance
#' @references M. G. B. Blum and O. Francois (2006). Which random processes
#' describe the Tree of Life? A large-scale study of phylogenetic tree
#' imbalance. Systematic Biology. 55:685-691.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#'   balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#'   unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#'   blum(balanced_tree)
#'   blum(unbalanced_tree) # should be higher
blum <- function(phy,
                 normalization = "none") {
  normalization <- check_normalization_key(normalization)
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    return(calc_blum_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_blum_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}


#' s shape statistic of (im)balance.
#' @description The s shape statistic of imbalance (also known as the Blum
#' statistic, see \link{blum}) calculates the sum of \eqn{log(N-1)} over all
#' internal nodes, where N represents the total number of extant tips connected
#' to that node. An alternative implementation can be found in the Castor R
#' package.
#' @param phy phylogeny or ltable
#' @param normalization "none" or "yule", in which case the acquired result
#' is divided by the expectation for the Yule model.
#' @return s shape statistic of imbalance
#' @references M. G. B. Blum and O. Francois (2006). Which random processes
#' describe the Tree of Life? A large-scale study of phylogenetic tree
#' imbalance. Systematic Biology. 55:685-691.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#'   balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#'   unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#'   sshape(balanced_tree)
#'   sshape(unbalanced_tree) # should be higher
sshape <- function(phy,
                   normalization = "none") {
  return(blum(phy, normalization))
}
