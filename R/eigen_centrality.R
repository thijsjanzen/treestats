#' Eigen vector centrality
#' @description Eigen vector centrality associates with each node \eqn{v}
#' the positive value \eqn{e(v)}, such that:
#' \eqn{ \sum_{e}^v w(uv) * e(u) = \lambda * e(v) }. Thus,
#' \eqn{e(v)} is the Perron-Frobenius eigenvector of the adjacency matrix of the
#' tree.
#' @param phy phylo object or ltable
#' @param weight if TRUE, uses branch lengths.
#' @param scale if TRUE, the Eigenvector is rescaled
#' @param use_rspectra boolean to indicate whether the helping package RSpectra
#' should be used, which is faster, but returns fewer eigen values.
#' @return List with the Eigen vector and the leading Eigen value
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
eigen_centrality <- function(phy,
                             weight = TRUE,
                             scale = FALSE,
                             use_rspectra = TRUE) {
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {

    if (use_rspectra == FALSE || !requireNamespace("RSpectra")) {
      ev <- slow_eigen_centrality(phy, weight)
    } else {
      edges   <- phy$edge
      lengths <- phy$edge.length
      if (!weight) {
        lengths <- rep(1, length(phy$edge.length))
      }
      n_nodes <- max(edges)

      matvec_fun <- function(x, args) {
        Ax_tree(edges, lengths, x, n_nodes)
      }

      mat_size <- max(phy$edge)
      ev <- RSpectra::eigs_sym(matvec_fun,
                               n = mat_size,
                               k = 1,
                               which = "LM",
                               opts = list(retvec = TRUE))
    }

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


slow_eigen_centrality <- function(phy, weight) {
  adj_matrix <- prep_adj_mat(as.vector(t(phy$edge)),
                             as.vector(phy$edge.length),
                             weight)
  ev <- eigen(adj_matrix, symmetric = TRUE)
  return(ev)
}
