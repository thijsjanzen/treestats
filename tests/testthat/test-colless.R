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


  # now, using ltable:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]
  focal_ltab <- treestats::phylo_to_l(focal_tree)

  colless <- treestats::colless(focal_ltab)
  colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree))
  testthat::expect_equal(colless, colless_check)

  colless <- treestats::colless(focal_ltab, normalization = "yule")
  colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                        norm = "yule")
  testthat::expect_equal(colless, colless_check, tol = 1e-5)

  colless <- treestats::colless(focal_ltab, normalization = "pda")
  colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                        norm = "pda")
  testthat::expect_equal(colless, colless_check, tol = 1e-5)


  ## with extinct lineages:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.2)[[1]]

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
})

test_that("abuse", {
  wrong_object <- list("a" = 5)
  testthat::expect_error(
    treestats::colless(wrong_object),
    "input object has to be phylo or ltable")
})
