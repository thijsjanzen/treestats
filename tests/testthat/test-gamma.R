context("gamma")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 120,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  gammast <- treestats::gamma_statistic(focal_tree)
  gammast_check <- ape::gammaStat(focal_tree)
  gammast_check2 <- castor::gamma_statistic(focal_tree)
  gammast_check3 <- phytools::ltt(focal_tree, plot = FALSE, gamma = TRUE)$gamma

  testthat::expect_equal(gammast_check, gammast_check2)
  testthat::expect_equal(gammast_check2, gammast_check3)

  testthat::expect_equal(gammast, gammast_check, tolerance = 1e-4)



  focal_trees <- TreeSim::sim.bd.taxa(n = 10, numbsim = 3, lambda = 1, mu = 0)

  gammast1 <- treestats::gamma_statistic(focal_trees)
  class(focal_trees) <- "multiPhylo"
  gammast2 <- treestats::gamma_statistic(focal_trees)
  testthat::expect_equal(gammast1, gammast2)
})
