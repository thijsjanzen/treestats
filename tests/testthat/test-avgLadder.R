context("avgLadder")

test_that("usage", {

  if (requireNamespace("phyloTop")) {
    standard_test(treestats::avg_ladder, phyloTop::avgLadder)
  }
})

test_that("polytomies", {
  for (focal_tree in poly_trees) {
    res <- try_stat(focal_tree, treestats::avg_ladder)
    testthat::expect_true(is.na(res))
  }
})

test_that("abuse", {
  test_abuse(treestats::avg_ladder)
})
