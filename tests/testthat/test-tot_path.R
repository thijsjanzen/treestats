context("tot_path")

test_that("usage", {

  if (requireNamespace("treebalance")) {
    standard_test(treestats::tot_path_length, treebalance::totPathLen)

    test_polytomies(treestats::tot_path_length, treebalance::totPathLen)
  }
})


test_that("yule", {
  testthat::skip_on_cran()
  test_yule(treestats::tot_path_length)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::tot_path_length(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::tot_path_length(list()),
    "input object has to be phylo or ltable"
  )
})
