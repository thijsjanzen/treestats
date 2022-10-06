context("avgLadder")

test_that("usage", {

  if (requireNamespace("phyloTop")) {
    set.seed(42)

    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    c1 <- treestats::avg_ladder(focal_tree)
    c2 <- phyloTop::avgLadder(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::avg_ladder(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)


    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0.5, fossils = TRUE)

    c1 <- treestats::avg_ladder(focal_tree)
    c2 <- phyloTop::avgLadder(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::avg_ladder(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)
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
