context("cherries")

test_that("usage", {
  if (requireNamespace("phyloTop")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    c1 <- treestats::cherries(focal_tree)
    c2 <- phyloTop::cherries(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::cherries(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)


    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0.5, fossils = TRUE)

    c1 <- treestats::cherries(focal_tree)
    c2 <- phyloTop::cherries(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::cherries(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)
  }
})

test_that("normalisation", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::cherries(focal_tree)
  c2 <- treestats::cherries(focal_tree, normalization = "yule")
  c3 <- treestats::cherries(focal_tree, normalization = "pda")
  testthat::expect_lt(c2, c1)
  testthat::expect_lt(c3, c1)

  focal_ltab <- treestats::phylo_to_l(focal_tree)

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
