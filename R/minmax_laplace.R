#' Laplacian Matrix properties
#' @description Calculates the eigenvalues of the Laplacian Matrix, where the
#' Laplacian matrix is the matrix representation of a graph, in this case a
#' phylogeny.
#' When the R package \pkg{RSpectra} is available, a faster
#' calculation can be used, which does not calculate all eigenvalues, but only
#' the maximum and minimum. As such, when using this option, the vector of all
#' eigenvalues is not returned.
#' @param phy phylo object or ltable
#' @param use_rspectra boolean to indicate whether the helping package RSpectra
#' should be used, in which case only the minimum and maximum values are
#' returned
#' @param force_min boolean to indicate whether or not to also calculate
#' the minimum positive eigenvalue if the tree is larger than 23170 tips.
#' Default is FALSE, as this can be very slow, but in some cases might be
#' of interest.
#' @return List with the minimum and maximum eigenvalues
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
minmax_laplace <- function(phy,
                           use_rspectra = FALSE,
                           force_min = FALSE) {
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  if (inherits(phy, "matrix")) {
    if (sum(phy[, 4] > -1))
      warning("minmax_laplace may give incorrect results for the minimal value due to conversion to a phylo object, when extinct lineages are present")  # nolint
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {

    if (!use_rspectra || !requireNamespace("RSpectra")) {
      return(slow_calc_minmax_laplace(phy))
    }

    edges   <- phy$edge
    lengths <- phy$edge.length
    n_nodes <- max(edges)

    matvec_fun <- function(x, args) {
      Lx_tree_weighted(edges, lengths, x, n_nodes)
    }

    mat_size <- max(phy$edge)
    max_val <- RSpectra::eigs_sym(matvec_fun,
                                  n = mat_size,
                                  k = 1,
                                  which = "LA",
                                  opts = list(retvec = FALSE))$values

    if (mat_size > 46340) {  # floor(sqrt(2^31 - 1)))
      if (force_min == TRUE) {
        warning("calculating the minimum eigenvalue
                   is very costly for large trees")
        res <- RSpectra::eigs_sym(matvec_fun,
                                  n = n_nodes,
                                  # we don't need negative eigenvalues:
                                  k = n_nodes / 2,
                                  which = "SM",
                                  opts = list(retvec = FALSE))$values
        res <- res[round(res, 10) > 0]
        res_min <- min(res)
        return(list("max" = max_val,
                    "min" = res_min))
      } else {
        warning("skipped calculating the minimum eigenvalue of the adjacency
                  matrix, because the tree is too large")
        return(list("max" = max_val))
      }
    } else {
      # if the tree is not very large, then using the matrix is faster
      lap_mat <- matrix(0, mat_size, mat_size)
      lap_mat[phy$edge] <- phy$edge.length
      lap_mat <- lap_mat + t(lap_mat)

      # we can make use of RSpectra
      min_val <- RSpectra::eigs_sym(lap_mat, k = 1, which = "BE",
                                    sigma = 1e-9,
                                    opts = list(retvec = FALSE))$values

      return(list(min = min_val,
                  max = max_val))
    }
  }
  stop("input object has to be phylo or ltable")
}




#' @keywords internal
slow_calc_minmax_laplace <- function(phy) {

  mat_size <- max(phy$edge)
  if (mat_size > 46340) {  # floor(sqrt(2^31 - 1)))
    stop("tree too big, memory allocation fail")
  }

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

