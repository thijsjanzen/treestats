context("mw_over_md")

test_that("usage", {

  if (requireNamespace("treebalance")) {
    standard_test(treestats::mw_over_md, treebalance::mWovermD)
    test_polytomies(treestats::mw_over_md, treebalance::mWovermD)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::mw_over_md(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mw_over_md(list()),
    "input object has to be phylo or ltable"
  )
})
