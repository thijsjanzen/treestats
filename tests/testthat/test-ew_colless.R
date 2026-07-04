context("ew_colless")

test_that("usage", {
  if (requireNamespace("treebalance")) {

    colless <- treestats::ew_colless(extant_tree)
    colless_check <- treebalance::ewCollessI(extant_tree)
    testthat::expect_equal(colless, colless_check)

    focal_ltab <- treestats::phylo_to_l(extant_tree)

    colless <- treestats::ew_colless(focal_ltab)
    testthat::expect_equal(colless, colless_check)

    ## with extinct lineages:
    colless <- treestats::ew_colless(complete_tree)
    colless_check <- treebalance::ewCollessI(complete_tree)
    testthat::expect_equal(colless, colless_check)

    focal_ltab <- treestats::phylo_to_l(complete_tree)

    colless <- treestats::ew_colless(focal_ltab)
    testthat::expect_equal(colless, colless_check)
  }
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {
    test_polytomies(treestats::ew_colless, treebalance::ewCollessI)
  }
})


test_that("wrong_object", {
  testthat::expect_error(
    treestats::ew_colless(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::ew_colless(list()),
    "input object has to be phylo or ltable"
  )
})
