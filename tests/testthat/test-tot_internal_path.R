context("tot_internal_path")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::tot_internal_path, treebalance::totIntPathLen)

    test_polytomies(treestats::tot_internal_path, treebalance::totIntPathLen)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::tot_internal_path(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::tot_internal_path(list()),
    "input object has to be phylo or ltable"
  )
})
