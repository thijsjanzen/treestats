#' Adjancency Matrix properties
#' @description Calculates the eigenvalues of the Adjancency Matrix, where the
#' Adjacency matrix is a square matrix indicate whether pairs of vertices
#' are adjacent or not on a graph - here, entries in the matrix indicate
#' connections between nodes (and betweens nodes and tips). Entries in the
#' adjacency matrix are weighted by branch length. Then, using the adjacency
#' matrix, we calculate the spectral properties of the matrix, e.g. the
#' minimum and maximum eigenvalues of the matrix.
#' @param phy phylo object or ltable
#' @return List with the minimum and maximum eigenvalues
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
minmax_adj <- function(phy) {
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    # prepping matrix in Rcpp yields no speed gain

    mat_size <- max(phy$edge)
    if (mat_size > 46340) {  # floor(sqrt(2^31 - 1)))
      Rcpp::stop("tree too big, memory allocation fail")
    }

    lap_mat <- matrix(0, mat_size, mat_size)
    lap_mat[phy$edge] <- ifelse(rep(TRUE, mat_size - 1), # weight is always true
                                phy$edge.length, 1)

    lap_mat <- lap_mat + t(lap_mat)

    eigen_vals <- eigen(lap_mat,
                        symmetric = TRUE,
                        only.values = TRUE)$values

    eigen_vals <- round(eigen_vals, digits = 10)

    max_val <- max(eigen_vals)
    min_val <- min(eigen_vals[eigen_vals > 0])
    return(list("max" = max_val,
                "min" = min_val,
                "values" = eigen_vals))
  }

  stop("input object has to be phylo or ltable")
}
