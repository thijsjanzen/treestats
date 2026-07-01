context("avgLadder")

test_that("usage", {

  if (requireNamespace("phyloTop")) {
    c1 <- treestats::avg_ladder(extant_tree)
    c2 <- phyloTop::avgLadder(extant_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::avg_ladder(treestats::phylo_to_l(extant_tree))
    testthat::expect_equal(c1, c3)

    c1 <- treestats::avg_ladder(complete_tree)
    c2 <- phyloTop::avgLadder(complete_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::avg_ladder(treestats::phylo_to_l(complete_tree))
    testthat::expect_equal(c1, c3)
  }
})

test_that("polytomies", {
  if (requireNamespace("phyloTop")) {
    test_polytomies(treestats::avg_ladder, phyloTop::avgLadder)
  }
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
