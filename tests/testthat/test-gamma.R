context("gamma")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 10, numbsim = 1, lambda = 1, mu = 0)[[1]]

  gammast <- treestats::gamma_statistic(focal_tree)
  gammast_check <- ape::gammaStat(focal_tree)
  testthat::expect_equal(gammast, gammast_check)
})
