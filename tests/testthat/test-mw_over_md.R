context("mw_over_md")

test_that("usage", {

  if (requireNamespace("treebalance")) {
    standard_test(treestats::mw_over_md, treebalance::mWovermD)
    test_polytomies(treestats::mw_over_md, treebalance::mWovermD)
  }
})

test_that("wrong_object", {
  test_abuse(treestats::mw_over_md)
})

