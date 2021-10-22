context("gamma")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]]

  gamma <- treestats::gamma(focal_tree)
  gamma_check <- ape::gammaStat(focal_tree)
  testthat::expect_equal(gamma, gamma_check)
})
