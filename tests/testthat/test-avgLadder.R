context("avgLadder")

test_that("usage", {

  if (requireNamespace("phyloTop")) {
    set.seed(42)

    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    c1 <- treestats::avgLadder(focal_tree)
    c2 <- phyloTop::avgLadder(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::avgLadder(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)


    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0.5, fossils = TRUE)

    c1 <- treestats::avgLadder(focal_tree)
    c2 <- phyloTop::avgLadder(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::avgLadder(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)
  }
})
