context("tot_internal_path")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::tot_internal_path, treebalance::totIntPathLen)

    test_polytomies(treestats::tot_internal_path, treebalance::totIntPathLen)
  }
})

test_that("yule", {
  testthat::skip_on_cran()
  test_yule(treestats::tot_internal_path)
})

test_that("wrong_object", {
  test_abuse(treestats::tot_internal_path)
})
