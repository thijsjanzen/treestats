#' Calculate the average leaf depth statistic. The average leaf depth statistic
#' is a normalized version of the Sackin index, normalized by the number of
#' tips.
#' @param phy phylo object or ltable
#' @return gamma statistic
#' @export
#' @references K.-T. Shao and R. R. Sokal. Tree balance. Systematic Zoology,
#' 39(3):266, 1990. doi: 10.2307/2992186.
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' average_leaf_depth(simulated_tree)
average_leaf_depth <- function(phy) {

  if (inherits(phy, "phylo")) {
    n <- length(phy$tip.label)
    return(calc_sackin_cpp(as.vector(t(phy$edge)),
                           normalization = "none") / n)
  }

  if (inherits(phy, "matrix")) {
    n <- length(phy[, 1])
    return(calc_sackin_ltable_cpp(phy,
                                  normalization = "none") / n)
  }

  stop("input object has to be phylo or ltable")
}
