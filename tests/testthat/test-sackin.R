context("sackin")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]]

  sackin <- treestats::sackin(focal_tree)
  sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(focal_tree))
  testthat::expect_equal(sackin, sackin_check)
})
