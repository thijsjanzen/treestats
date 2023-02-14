context("max_ladder")

test_that("usage", {

  if (requireNamespace("phyloTop")) {
    set.seed(42)

    focal_tree <- ape::rphylo(n = 1000, birth = 1, death = 0)

    c1 <- treestats::max_ladder(focal_tree)
    c2 <- max(phyloTop::ladderSizes(focal_tree)$ladderSizes)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::max_ladder(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)

    # we use larger trees here to induce larger ladders.
    focal_tree <- ape::rphylo(n = 1000, birth = 1, death = 0.5, fossils = TRUE)

    c1 <- treestats::max_ladder(focal_tree)
    c2 <- max(phyloTop::ladderSizes(focal_tree)$ladderSizes)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::max_ladder(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::max_ladder(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::max_ladder(list()),
    "input object has to be phylo or ltable"
  )
})
