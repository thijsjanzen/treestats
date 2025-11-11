#' Colless index of (im)balance.
#' @description The Colless index is calculated as the sum of
#' \eqn{abs(L - R)} over all nodes, where L (or R) is the number of extant tips
#' associated with the L (or R) daughter branch at that node.  Higher values
#' indicate higher imbalance. Two normalizations are available,
#' where a correction is made for tree size, under either a yule expectation,
#' or a pda expectation.
#' @param phy phylo object or ltable
#' @param normalization A character string equals to "none" (default) for no
#' normalization or one of "pda" or "yule".
#' @return colless index
#' @references  Colless D H. 1982. Review of: Phylogenetics: The Theory and
#' Practice of Phylogenetic Systematics. Systematic Zoology 31:100-104.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#' unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#' colless(balanced_tree)
#' colless(unbalanced_tree) # should be higher
colless <- function(phy,
                    normalization = "none") {
  normalization <- check_normalization_key(normalization)
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    return(calc_colless_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_colless_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}

#' Equal weights Colless index of (im)balance.
#' @description The equal weights Colless index is calculated as the sum of
#' \eqn{abs(L - R) / (L + R - 2)} over all nodes where \eqn{L + R > 2},
#' where L (or R) is the number of extant tips associated with the L (or R)
#' daughter branch at that node.  Maximal imbalance is associated with a value
#' of 1.0. The ew_colless index is not sensitive to tree size.
#' @param phy phylo object or ltable
#' @return colless index
#' @references  A. O. Mooers and S. B. Heard. Inferring Evolutionary Process
#' from Phylogenetic Tree Shape. The Quarterly Review of Biology, 72(1), 1997.
#' doi: 10.1086/419657.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#' unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#' ew_colless(balanced_tree)
#' ew_colless(unbalanced_tree) # should be higher
ew_colless <- function(phy) {
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)
  if (inherits(phy, "matrix")) {
    return(calc_eWcolless_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_eWcolless_cpp(as.vector(t(phy$edge))))
  }
  stop("input object has to be phylo or ltable")
}



#' Corrected Colless index of (im)balance.
#' @description The Corrected Colless index is calculated as the sum of
#' \eqn{abs(L - R)} over all nodes, corrected for tree size by dividing over
#' \eqn{(n-1) * (n-2)}, where n is the number of nodes.
#' @param phy phylo object or ltable
#' @param normalization A character string equals to "none" (default) for no
#' normalization or "yule", in which case the obtained index is divided by
#' the Yule expectation.
#' @return corrected colless index
#' @references  Heard, Stephen B. "Patterns in tree balance among cladistic,
#' phenetic, and randomly generated phylogenetic trees." Evolution 46.6 (1992):
#' 1818-1826.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#' unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#' colless_corr(balanced_tree)
#' colless_corr(unbalanced_tree) # should be higher
colless_corr <- function(phy,
                         normalization = "none") {
  normalization <- check_normalization_key(normalization)
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    return(calc_colless_corr_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_colless_corr_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}

#' Quadratic Colless index of (im)balance.
#' @description The Quadratic Colless index is calculated as the sum of
#' \eqn{(L - R)^2} over all nodes.
#' @param phy phylo object or ltable
#' @param normalization A character string equals to "none" (default) for no
#' normalization or "yule"
#' @return quadratic colless index
#' @references  Bartoszek, Krzysztof, et al. "Squaring within the Colless index
#' yields a better balance index." Mathematical Biosciences 331 (2021): 108503.
#' @export
#' @examples simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
#' balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
#' unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
#' colless_quad(balanced_tree)
#' colless_quad(unbalanced_tree) # should be higher
colless_quad <- function(phy,
                         normalization = "none") {
  normalization <- check_normalization_key(normalization)
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    return(calc_colless_quad_ltable_cpp(phy, normalization))
  }
  if (inherits(phy, "phylo")) {
    return(calc_colless_quad_cpp(as.vector(t(phy$edge)), normalization))
  }
  stop("input object has to be phylo or ltable")
}
