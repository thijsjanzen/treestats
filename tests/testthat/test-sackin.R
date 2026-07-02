context("sackin")

test_that("usage", {
  if (requireNamespace("treebalance")) {

    standard_test(treestats::sackin, treebalance::sackinI)

    test_polytomies(treestats::sackin, treebalance::sackinI)
  }

  if (requireNamespace("apTreeshape")) {

    test_func <- function(phy) {
      apTreeshape::sackin(apTreeshape::as.treeshape(phy))
    }

    standard_test(treestats::sackin, test_func)

    # no polytomies test, as apTreeshape does not support it.

    sackin <- treestats::sackin(extant_tree, normalization = "yule")
    sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(extant_tree),
                                        norm = "yule")
    testthat::expect_equal(sackin, sackin_check, tol = 1e-5)

    sackin <- treestats::sackin(extant_tree, normalization = "pda")
    sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(extant_tree),
                                        norm = "pda")
    testthat::expect_equal(sackin, sackin_check, tol = 1e-5)


    # test ltable functionality
    ltab <- treestats::phylo_to_l(extant_tree)
    testthat::expect_equal(treestats::sackin(extant_tree),
                           treestats::sackin(ltab))

    testthat::expect_equal(treestats::sackin(extant_tree,
                                             normalization = "yule"),
                           treestats::sackin(ltab,
                                             normalization = "yule"))

    testthat::expect_equal(treestats::sackin(extant_tree,
                                             normalization = "pda"),
                           treestats::sackin(ltab,
                                             normalization = "pda"))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::sackin(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::sackin(list()),
    "input object has to be phylo or ltable"
  )
})
