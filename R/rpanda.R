#' @keywords internal
integr <- function(x, f) {
  # get lengths of var and integrand
  n <- length(x)

  # trapezoidal integration
  integral <- 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] + f[1:(n - 1)]))

  # print definite integral
  return(integral)
}


#' @keywords internal
skew_ness <- function(x) {
  n <- length(x)
  return((sum((x - mean(x)) ^ 3) / n) / (sum((x - mean(x)) ^ 2) / n) ^ (3 / 2))
}

#' Laplacian spectrum statistics, from RPANDA
#' @description Computes the distribution of eigenvalues for the modified graph
#' Laplacian of a phylogenetic tree, and several summary statistics of this
#' distribution. The modified graph Laplacian of a phylogeny is given by the
#' difference between its' distance matrix (e.g. all pairwise distances between
#' all nodes), and the degree matrix (e.g. the diagonal matrix where each
#' diagonal element represents the sum of branch lengths to all other nodes).
#' Each row of the modified graph Laplacian sums to zero. For a tree with n
#' tips, there are \eqn{N = 2n-1} nodes, and hence the modified graph Laplacian
#' is represented by a \eqn{N x N} matrix. Where RPANDA relies on the package
#' igraph to calculate the modified graph Laplacian, the treestats package uses
#' C++ to directly calculate the different entries in the matrix. This makes the
#' treestats implementation slightly faster, although the bulk of computation
#' occurs in estimating the eigen values, where RcppArmadillo generates a bit
#' of speed up.
#' @param phy phy
#' @return list with five components: 1) eigenvalues the vector of eigen
#' values, 2) principal_eigenvalue the largest eigenvalueof the spectral
#' density distribution 3) asymmetry the skewness of the spectral density
#' distribution 4) peak_height the largest y-axis valueof the spectral
#' density distribution and 5) eigengap theposition ofthe largest
#' difference between eigenvalues, giving the number of modalities in the tree.
#' @references Eric Lewitus, Helene Morlon, Characterizing and Comparing
#' Phylogenies from their Laplacian Spectrum, Systematic Biology, Volume 65,
#' Issue 3, May 2016, Pages 495–507, https://doi.org/10.1093/sysbio/syv116
#' @export
laplacian_spectrum <- function(phy) {

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (!inherits(phy, "phylo")) {
    stop("input object has to be phylo or ltable")
  }
  n <- length(phy$tip.label)
  m <- n - 1
  nm <- n + m
  if (nm > 46340) { # sqrt(2^31 - 1) #nolint
    stop("tree too big")
  }

  dens_rpanda <- function(x,
                          bw = stats::bw.nrd0,
                          n = 4096,
                          from = min(x) - 3 * sd, to = max(x) + 3 * sd,
                          adjust = 1,
                          ...) {
    if (has_na <- any(is.na(x))) {
      x <- stats::na.omit(x)
      if (length(x) == 0)
        stop("no finite or non-missing data!")
    }
    sd <- (if (is.numeric(bw))
      bw[1]
      else bw(x)) * adjust
    xx <- seq(from, to, len = n)

    mat <-  outer_cpp(xx, x, sd)


    structure(list(x = xx,
                   y = rowMeans(mat),
                   bw = sd,
                   call = match.call(),
                   n = length(x),
                   data.name = deparse(substitute(x)),
                   has.na = has_na),
              class = "density")
  }

  phy <- ape::reorder.phylo(phy)

  x <- get_eigen_values_arma_cpp(phy)
  x <- sort(x, decreasing = TRUE)

  d <- dens_rpanda(log(x))
  dsc <- d$y / (integr(d$x, d$y))
  principal_eigenvalue <- max(x)
  skewness <- skew_ness(x)
  peak_height <- max(dsc)
  gaps <- abs(diff(x))
  gap_mat <- as.matrix(gaps)  # nolint
  modalities <- seq_along(gap_mat)
  gap_mat_col <- cbind(modalities, gap_mat)  # nolint
  eigen_gap <- subset(gap_mat_col, gap_mat_col[, 2] == max(gap_mat_col[, 2]))
  res <- list(eigenvalues = x,
              principal_eigenvalue = principal_eigenvalue,
              asymmetry = skewness,
              peakedness = peak_height,
              eigengap = eigen_gap[, 1])
  return(res)
}
