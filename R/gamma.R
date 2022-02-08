#' Calculate the gamma statistic, using a fast implementation in C++.
#' @description The gamma statistic measures the relative position of
#' internal nodes within a reconstructed phylogeny. Under the Yule process,
#' the gamma values of a reconstructed tree follow a standard normal
#' distribution. If gamma > 0, the nodes are located more towards the tips of
#' the tree, and if gamma < 0, the nodes are located more towards the root of
#' the tree.
#' @param phy phylo object
#' @return gamma statistic
#' @export
#' @references Pybus, O. G. and Harvey, P. H. (2000) Testing macro-evolutionary
#' models using incomplete molecular phylogenies. Proceedings of the Royal
#' Society of London. Series B. Biological Sciences, 267, 2267–2272.
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' gamma_statistic(simulated_tree) # should be around 0.
#' ddd_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 7)$tes
#' gamma_statistic(ddd_tree) # because of diversity dependence, should be < 0
gamma_statistic <- function(phy) {

  apply_function_phy_ltable(phy,
                            calc_gamma_cpp,
                            calc_gamma_ltable_cpp,
                            only_extant = TRUE)
}
