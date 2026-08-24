context("stairs")

test_that("usage stairs1", {
  if (requireNamespace("phyloTop")) {
    test_func <- function(phy) {
      phyloTop::stairs(phy)[[1]]
    }

    standard_test(treestats::stairs, test_func)
  }
})

test_that("usage stairs2", {
  if (requireNamespace("phyloTop")) {

    test_func <- function(phy) {
      phyloTop::stairs(phy)[[2]]    # second index is stairs2
    }

    standard_test(treestats::stairs2, test_func)
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

test_that("polytomies", {
  for (focal_tree in poly_trees) {
    a1 <- try_stat(focal_tree, treestats::stairs)
    testthat::expect_true(is.na(a1))
    a1 <- try_stat(focal_tree, treestats::stairs2)
    testthat::expect_true(is.na(a1))
  }
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
