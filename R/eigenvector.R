#' Fast function using C++ to calculate eigen vector centrality
#' @description Eigen vector centrality associates with each node v the positive
#' value e(v), such that: \eqn{sum_{e~v} w(uv) * e(u) = \lambda * e(v)}. Thus,
#' e(v) is the Perron-Frobenius eigenvector of the adjacency matrix of the tree.
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @param scale if TRUE, the eigenvector is rescaled
#' @return List with the eigen vector and the leading eigen value
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
eigen_vector <- function(phy,
                         weight = TRUE,
                         scale = FALSE) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    edge_for_mat <- rbind(phy$edge, cbind(phy$edge[, 2], phy$edge[, 1]))

    adj_matrix <- c()
    if (weight) {
      adj_matrix <- Matrix::sparseMatrix(i = edge_for_mat[, 1],
                                         j = edge_for_mat[, 2],
                                         x = c(phy$edge.length,
                                                   phy$edge.length))
    } else {
      adj_matrix <- Matrix::sparseMatrix(i = edge_for_mat[, 1],
                                       j = edge_for_mat[, 2],
                                       x = rep(1, length(edge_for_mat[, 1])))
    }

    ev <- RSpectra::eigs_sym(adj_matrix, k = 1,
                             which = "LM",
                             opts = list(retvec = TRUE))
    evector <- abs(ev$vectors[, 1])

    evalue <- abs(ev$values[1])

    if (scale) {
      evector <- evector / max(evector)
    }

    return(list(eigenvector = evector,
                eigenvalue = evalue))
  }

  stop("input object has to be phylo or ltable")
}
