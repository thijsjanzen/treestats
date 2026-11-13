context("b2")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::b2, treebalance::B2I)
  }
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {
    test_polytomies(treestats::b2, treebalance::B2I)
  }
})

test_that("normalisation", {
  test_yule(treestats::b2)
})

test_that("abuse", {
 test_abuse(treestats::b2)
})
