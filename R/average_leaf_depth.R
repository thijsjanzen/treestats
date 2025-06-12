#' Average leaf depth statistic. The average leaf depth statistic
#' is a normalized version of the Sackin index, normalized by the number of
#' tips.
#' @param phy phylo object or ltable
#' @param normalization "none" or "yule", in which case the statistic is
#' divided by the expectation under the yule model, following Remark 1 in
#' Coronado et al. 2020.
#' @return average leaf depth statistic
#' @export
#' @references M. Coronado, T., Mir, A., Rosselló, F. et al. On Sackin’s
#' original proposal: the variance of the leaves’ depths as a phylogenetic
#' balance index. BMC Bioinformatics 21, 154 (2020).
#' https://doi.org/10.1186/s12859-020-3405-1
#' K.-T. Shao and R. R. Sokal. Tree balance. Systematic Zoology,
#' 39(3):266, 1990. doi: 10.2307/2992186.
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' average_leaf_depth(simulated_tree)
average_leaf_depth <- function(phy, normalization = "none") {
  normalization <- check_normalization_key(normalization)
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)


  if (inherits(phy, "phylo")) {
    n <- length(phy$tip.label)
    ald <- calc_sackin_cpp(as.vector(t(phy$edge)), "none") / n
    if (normalization == "yule" || normalization == TRUE) {
      expectation <- 2 * log(n)
      ald <- ald / expectation
    }
    return(ald)
  }

  if (inherits(phy, "matrix")) {
    n <- length(phy[, 1])
    ald <- calc_sackin_ltable_cpp(phy, "none") / n
    if (normalization == "yule" || normalization == TRUE) {
      expectation <- 2 * log(n)
      ald <- ald / expectation
    }
    return(ald)
  }

  stop("input object has to be phylo or ltable")
}
