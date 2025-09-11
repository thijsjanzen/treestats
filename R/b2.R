#' B2 metric
#' @description Balance metric that uses the Shannon-Wiener statistic of
#' information content. The b2 measure is given by the sum over the depths of
#' all tips, divided by \eqn{2^{d}}, e.g. \eqn{\sum d_i / 2^{d_i}}
#' Although also defined on non-binary trees, the treestats package only
#' provides code for binary trees.
#' @param phy phylo object or ltable
#' @param normalization "none" or "yule", when "yule" is chosen, the statistic
#' is divided by the Yule expectation, following from theorem 3.7 in Bienvenu
#' 2020.
#' @return Maximum depth (in number of edges)
#' @references  K.-T. Shao and R. R. Sokal. Tree Balance.
#' Systematic Zoology, 39(3):266, 1990. doi: 10.2307/2992186.
#'
#' Bienvenu, François, Gabriel Cardona, and Celine Scornavacca.
#' "Revisiting Shao and Sokal’s $$ B_2 $$ B 2 index of phylogenetic balance."
#' Journal of Mathematical Biology 83.5 (2021): 1-43.
#' @export
b2 <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    b2_stat <- calc_b2_ltable_cpp(phy)
    if (normalization == "yule" || normalization == TRUE) {
      n <- length(phy[, 1])
      expectation <- sum(1 / (1:n))
      b2_stat <- b2_stat / expectation
    }
    return(b2_stat)
  }
  if (inherits(phy, "phylo")) {
    b2_stat <- calc_b2_cpp(as.vector(t(phy$edge)))
    if (normalization == "yule" || normalization == TRUE) {
      n <- length(phy$tip.label)
      expectation <- sum(1 / (1:n))
      b2_stat <- b2_stat / expectation
    }
    return(b2_stat)
  }
  stop("input object has to be phylo or ltable")
}
