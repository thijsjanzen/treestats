context("pigot_rho")

test_that("usage", {
  set.seed(42)

  # DDD tree expected to slow down diversification
  focal_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 7)$tes
  rho <- treestats::pigot_rho(focal_tree)
  testthat::expect_lt(rho, 0)

  set.seed(42)
  # BD tree, expected increase in diversification
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.0)[[1]]
  rho <- treestats::pigot_rho(focal_tree)
  testthat::expect_gt(rho, 0)
})
