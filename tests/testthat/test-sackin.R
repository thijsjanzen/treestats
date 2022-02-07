context("sackin")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  sackin <- treestats::sackin(focal_tree)
  sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(focal_tree))
  testthat::expect_equal(sackin, sackin_check)

  sackin <- treestats::sackin(focal_tree, normalization = "yule")
  sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(focal_tree),
                                      norm = "yule")
  testthat::expect_equal(sackin, sackin_check, tol = 1e-5)

  sackin <- treestats::sackin(focal_tree, normalization = "pda")
  sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(focal_tree),
                                      norm = "pda")
  testthat::expect_equal(sackin, sackin_check, tol = 1e-5)
})

test_that("abuse", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.3)[[1]]
  testthat::expect_error(
    treestats::sackin(focal_tree),
    "can only calculate sackin statistic for ultrametric tree"
  )
})
