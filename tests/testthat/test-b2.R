context("b2")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::b2, treebalance::B2I)
  }
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {
    # TODO: this seems to crash for no reason?
    # test_polytomies(treestats::b2, treebalance::b2I)
  }
})

test_that("normalisation", {
  c1 <- treestats::b2(extant_tree)
  c2 <- treestats::b2(extant_tree, normalization = "yule")
  testthat::expect_lt(c2, c1)
  c3 <- treestats::b2(treestats::phylo_to_l(extant_tree),
                                    normalization = "yule")
  testthat::expect_equal(c2, c3)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::b2(focal_tree))
    stats2 <- c(stats2, treestats::b2(focal_tree, normalization = "yule"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.5)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::b2(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::b2(list()),
    "input object has to be phylo or ltable"
  )
})
