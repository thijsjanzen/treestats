#' Gamma statistic
#' @description The gamma statistic measures the relative position of
#' internal nodes within a reconstructed phylogeny. Under the Yule process,
#' the gamma values of a reconstructed tree follow a standard normal
#' distribution. If gamma > 0, the nodes are located more towards the tips of
#' the tree, and if gamma < 0, the nodes are located more towards the root of
#' the tree. Only available for ultrametric trees.
#' @param phy phylo object or ltable
#' @return gamma statistic
#' @export
#' @references Pybus, O. G. and Harvey, P. H. (2000) Testing macro-evolutionary
#' models using incomplete molecular phylogenies. Proceedings of the Royal
#' Society of London. Series B. Biological Sciences, 267, 2267–2272.
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' gamma_statistic(simulated_tree) # should be around 0.
#' if (requireNamespace("DDD")) {
#'   ddd_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 7)$tes
#'   gamma_statistic(ddd_tree) # because of diversity dependence, should be < 0
#' }
gamma_statistic <- function(phy) {

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = TRUE,
             require_rooted = TRUE)

  if (inherits(phy, "phylo")) {

    # if (length(phy$tip.label) < 100) {
    #  return(calc_gamma_cpp(phy))
    # } else {
    #   return(calc_gamma_cpp2(as.vector(t(phy$edge)), phy$edge.length))
    # }
    return(calc_gamma_cpp(phy))
  }
  if (inherits(phy, "matrix")) {
    return(calc_gamma_ltable_cpp(phy))
  }

  stop("input object has to be phylo or ltable")
}
