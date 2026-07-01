context("stairs")

test_that("usage stairs1", {
  if (requireNamespace("phyloTop")) {

    test_func <- function(phy) {
      phyloTop::stairs(phy)[[1]]
    }

    standard_test(treestats::stairs, test_func)

    test_polytomies(treestats::stairs, test_func)
  }
})

test_that("usage stairs2", {
  if (requireNamespace("phyloTop")) {

    test_func <- function(phy) {
      phyloTop::stairs(phy)[[2]] # 2 = stairs2
    }

    standard_test(treestats::stairs2, test_func)

    test_polytomies(treestats::stairs2, test_func)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::stairs(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::stairs(list()),
    "input object has to be phylo or ltable"
  )
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::stairs2(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::stairs2(list()),
    "input object has to be phylo or ltable"
  )
})
