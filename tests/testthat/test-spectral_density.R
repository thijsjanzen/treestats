context("spectral density")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  spect_dens <- treestats::spectral_density(focal_tree)
  spect_dens_check <- RPANDA::spectR(focal_tree)
  testthat::expect_equal(spect_dens, spect_dens_check)
})
