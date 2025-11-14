#' Adjancency Matrix properties
#' @description Calculates the eigenvalues of the Adjancency Matrix, where the
#' Adjacency matrix is a square matrix indicate whether pairs of vertices
#' are adjacent or not on a graph - here, entries in the matrix indicate
#' connections between nodes (and betweens nodes and tips). Entries in the
#' adjacency matrix are weighted by branch length. Then, using the adjacency
#' matrix, we calculate the spectral properties of the matrix, e.g. the
#' minimum positive and maximum positive eigenvalues of the matrix.
#' When the R package \pkg{RSpectra} is available,
#' a faster calculation can be used, which does not calculate all eigenvalues,
#' but only the maximum and minimum. As such, when using this option,
#' the vector of all eigenvalues is not returned.
#' @param phy phylo object or ltable
#' @param use_rspectra boolean to indicate whether the helping package RSpectra
#' should be used, in which case only the minimum and maximum values are
#' returned
#' @param calculate_min boolean to indicate whether or not to also calculate
#' the minimum positive eigenvalue. Default is TRUE, but when calculating for
#' extremely large trees (>10k tips, or even >40k tips), this can be very slow
#' or in which case setting to FALSE still allows for the
#' calculation of the maximum eigen value within reasonable time.
#' @return List with the minimum and maximum eigenvalues
#' @references  Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." Plos one 16.12 (2021): e0259877.
#' @export
minmax_adj <- function(phy, use_rspectra = FALSE, calculate_min = TRUE) {
  check_tree(phy,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = FALSE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {

    # when not using RSpectra, things are slow and restricted
    if (!use_rspectra || !requireNamespace("RSpectra")) {
      return(slow_calc(phy))
    } else {
      edges   <- phy$edge
      lengths <- phy$edge.length
      n_nodes <- max(edges)

      matvec_fun <- function(x, args) {
        Ax_tree(edges, lengths, x, n_nodes)
      }

      mat_size <- max(phy$edge)
      max_val <- RSpectra::eigs_sym(matvec_fun,
                                    n = mat_size,
                                    k = 1,
                                    which = "LA",
                                    opts = list(retvec = FALSE))$values

      if (mat_size > 46340) {  # floor(sqrt(2^31 - 1)))
        if (calculate_min == TRUE) {
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
  }

  stop("input object has to be phylo or ltable")
}


#' @keywords internal
slow_calc <- function(phy) {
  warning("calculating all eigenvalues without RSpectra can be slow
              for large trees")
  # prepping matrix in Rcpp yields no speed gain
  mat_size <- max(phy$edge)
  if (mat_size > 46340) {  # floor(sqrt(2^31 - 1)))
    stop("tree too big, memory allocation fail")
  }
  lap_mat <- matrix(0, mat_size, mat_size)
  lap_mat[phy$edge] <- phy$edge.length
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
