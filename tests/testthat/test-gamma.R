context("gamma")



test_that("usage", {

  check_gamma <- function(focal_tree) {
    gammast <- treestats::gamma_statistic(focal_tree)
    gammast_check <- ape::gammaStat(focal_tree)
    gammast_check2 <- castor::gamma_statistic(focal_tree)
    gammast_check3 <- phytools::ltt(focal_tree, plot = FALSE, gamma = TRUE)$gamma

    testthat::expect_equal(gammast, gammast_check, tolerance = 1e-4)
    testthat::expect_equal(gammast, gammast_check2, tolerance = 1e-4)
    testthat::expect_equal(gammast, gammast_check3, tolerance = 1e-4)
  }

  set.seed(42)
  bd_tree <- TreeSim::sim.bd.taxa(n = 121,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  check_gamma(bd_tree)

  # DDD tree has all speciation events at the beginning
  # expect strongly negative gamma.
  dd_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 20)$tas
  check_gamma(dd_tree)
  testthat::expect_lt(treestats::gamma_statistic(dd_tree), 0)
})
