context("colless")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  colless <- treestats::colless(focal_tree)
  colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree))
  testthat::expect_equal(colless, colless_check)

  colless <- treestats::colless(focal_tree, normalization = "yule")
  colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                      norm = "yule")
  testthat::expect_equal(colless, colless_check, tol = 1e-5)

  colless <- treestats::colless(focal_tree, normalization = "pda")
  colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                      norm = "pda")
  testthat::expect_equal(colless, colless_check, tol = 1e-5)

  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.2)[[1]]

  testthat::expect_error(
    treestats::colless(focal_tree)
  )

})
