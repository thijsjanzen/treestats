context("avgLadder")

test_that("usage", {

  if (requireNamespace("phyloTop")) {
    standard_test(treestats::avg_ladder, phyloTop::avgLadder)
  }
})

test_that("polytomies", {
  # phyloTop transforms the trees, e.g. forces them to binary,
  # which is a no-no.

  #if (requireNamespace("phyloTop")) {
  #  test_polytomies(treestats::avg_ladder, phyloTop::avgLadder)
  #}
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::avg_ladder(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::avg_ladder(list()),
    "input object has to be phylo or ltable"
  )
})
