context("rpanda")

test_that("usage", {
  skip_on_cran() # RPANDA uses iGraph for the Laplacian spectrum
                 # this can cause BLAS issues

  if (requireNamespace("RPANDA")) {
    set.seed(42)
    focal_tree <- TreeSim::sim.bd.taxa(n = 10,
                                       numbsim = 1,
                                       lambda = 1, mu = 0)[[1]]

    ref <- RPANDA::spectR(focal_tree)
    stat <- treestats::calc_lapl_spectrum(focal_tree)

    diff_eig <- sum(ref$eigenvalues - stat$eigenvalues)
    testthat::expect_equal(0, diff_eig)
    testthat::expect_equal(ref$asymmetry, stat$asymmetry)
    testthat::expect_equal(ref$peakedness, stat$peakedness, tolerance = 0.01)
    testthat::expect_equal(ref$principal_eigenvalue, stat$principal_eigenvalue)
    testthat::expect_equal(ref$eigengap, stat$eigengap)
  }
})
