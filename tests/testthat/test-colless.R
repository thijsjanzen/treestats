context("colless")

test_that("usage", {
  if (requireNamespace("apTreeshape")) {
     colless <- treestats::colless(extant_tree)
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(extant_tree))
    testthat::expect_equal(colless, colless_check)

    colless <- treestats::colless(extant_tree, normalization = "yule")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(extant_tree),
                                          norm = "yule")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)

    colless <- treestats::colless(extant_tree, normalization = "pda")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(extant_tree),
                                          norm = "pda")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)

    # now, using ltable:
    focal_ltab <- treestats::phylo_to_l(extant_tree)

    colless <- treestats::colless(focal_ltab)
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(extant_tree))
    testthat::expect_equal(colless, colless_check)

    colless <- treestats::colless(focal_ltab, normalization = "yule")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(extant_tree),
                                          norm = "yule")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)

    colless <- treestats::colless(focal_ltab, normalization = "pda")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(extant_tree),
                                          norm = "pda")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)


    ## with extinct lineages:
    colless <- treestats::colless(complete_tree)
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(complete_tree))
    testthat::expect_equal(colless, colless_check)

    colless <- treestats::colless(complete_tree, normalization = "yule")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(complete_tree),
                                          norm = "yule")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)

    colless <- treestats::colless(complete_tree, normalization = "pda")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(complete_tree),
                                          norm = "pda")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)
  }
})

test_that("polytomies", {
  if (requireNamespace("apTreeshape")) {

    test_func <- function(phy) {
      apTreeshape::colless(apTreeshape::as.treeshape(phy))
    }

    test_polytomies(treestats::colless, test_func)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::colless(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::colless(list()),
    "input object has to be phylo or ltable"
  )
})
