#' Laplacian Matrix properties
#' @description Calculates the eigenvalues of the Laplacian Matrix, where the
#' Laplacian matrix is the matrix representation of a graph, in this case a
#' phylogeny.
#' @param phy phylo object or ltable
#' @return List with the minimum and maximum eigenvalues
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
minmax_laplace <- function(phy) {
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(phy, "matrix")) {
    if (sum(phy[, 4] > -1))
    warning("minmax_laplace may give incorrect results for the minimal value due to conversion to a phylo object, when extinct lineages are present")  # nolint
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    mat_size <- max(phy$edge)
    lap_mat <- matrix(0, mat_size, mat_size)
    lap_mat[phy$edge] <- ifelse(rep(TRUE, mat_size - 1), # weight is always true
                                phy$edge.length, 1)

    lap_mat <- lap_mat + t(lap_mat)
    degrees <- rowSums(lap_mat)
    lap_mat <- diag(degrees) - lap_mat

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
