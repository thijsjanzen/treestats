#' Fast function using C++ to calculate the rquartet index.
#' @description The rquartet index counts the number of potential fully
#' balanced rooted subtrees of 4 tips in the tree. The function in treestats
#' assumes a bifurcating tree. For trees with polytomies, we refer the user to
#' treebalance::rquartedI, which can also take polytomies into account.
#' @param phy phylo object or ltable
#' @param normalization The index can be normalized by the expectation under
#' the Yule ("yule") or PDA model ("pda").
#' @return rquartet index
#' @references  T. M. Coronado, A. Mir, F. Rosselló, and G. Valiente. A balance
#' index for phylogenetic trees based on rooted quartets. Journal of
#' Mathematical Biology, 79(3):1105-1148, 2019. doi: 10.1007/s00285-019-01377-w.
#' @export
rquartet <- function(phy, normalization = "none") {
  if (inherits(phy, "matrix")) {
    answ <- calc_rquartet_ltable_cpp(phy)
    if (normalization == "yule" || normalization == "Yule") {
      answ <- answ / (choose(length(phy[, 1]), 4))
    }
    if (normalization == "pda" || normalization == "PDA") {
      answ <- answ / (5 * choose(length(phy[, 1]), 4))
    }
    return(answ)
  }
  if (inherits(phy, "phylo")) {
    if (phy$Nnode + 1 != length(phy$tip.label)) {
  stop("Tree must be binary, for non binary trees use treebalance::rQuartetI")
    }
    answ <- calc_rquartet_cpp(as.vector(t(phy$edge)))
    if (normalization == "yule" || normalization == "Yule") {
      answ <- answ / (choose(length(phy$tip.label), 4))
    }
    if (normalization == "pda" || normalization == "PDA") {
      answ <- answ / (5 * choose(length(phy$tip.label), 4))
    }
    return(answ)
  }
  stop("input object has to be phylo or ltable")
}
