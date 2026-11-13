context("var leaf depth")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::var_leaf_depth, treebalance::varLeafDepI)

    test_polytomies(treestats::var_leaf_depth, treebalance::varLeafDepI)
  }
})


test_that("normalization", {
  testthat::skip_on_cran()
  test_yule(treestats::var_leaf_depth)
})

test_that("wrong_object", {
  test_abuse(treestats::var_leaf_depth)
})
