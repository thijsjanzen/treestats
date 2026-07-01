context("cherries")

test_that("usage", {
  if (requireNamespace("phyloTop")) {
    c1 <- treestats::cherries(extant_tree)
    c2 <- phyloTop::cherries(extant_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::cherries(treestats::phylo_to_l(extant_tree))
    testthat::expect_equal(c1, c3)

    c1 <- treestats::cherries(complete_tree)
    c2 <- phyloTop::cherries(complete_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::cherries(treestats::phylo_to_l(complete_tree))
    testthat::expect_equal(c1, c3)
  }
})

test_that("normalisation", {
  c1 <- treestats::cherries(extant_tree)
  c2 <- treestats::cherries(extant_tree, normalization = "yule")
  c3 <- treestats::cherries(extant_tree, normalization = "pda")
  testthat::expect_lt(c2, c1)
  testthat::expect_lt(c3, c1)

  focal_ltab <- treestats::phylo_to_l(extant_tree)

  c4 <- treestats::cherries(focal_ltab)
  c5 <- treestats::cherries(focal_ltab, normalization = "yule")
  c6 <- treestats::cherries(focal_ltab, normalization = "pda")

  testthat::expect_equal(c1, c4)
  testthat::expect_equal(c2, c5)
  testthat::expect_equal(c3, c6)

  stats1 <- c()
  stats2 <- c()
  stats3 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::cherries(focal_tree))
    stats2 <- c(stats2, treestats::cherries(focal_tree,
                                            normalization = "yule"))
    stats3 <- c(stats3, treestats::cherries(focal_tree,
                                            normalization = "pda"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))
  a3 <- cor(stats3, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_equal(a2, a3)
})

test_that("polytomies", {
  # phyloTop does not support non-binary trees, so no comparison possible.
  for (focal_tree in poly_trees) {
   local_stat <- try_stat(focal_tree, treestats::cherries)
   testthat::expect_true(!is.na(local_stat))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::cherries(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::cherries(list()),
    "input object has to be phylo or ltable"
  )
})
