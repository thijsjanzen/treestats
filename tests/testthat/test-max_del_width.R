context("max del width")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::max_del_width(focal_tree)
    a2 <- treebalance::maxDelW(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::max_del_width(focal_tree),
                           treestats::max_del_width(ltab))

    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    a1 <- treestats::max_del_width(focal_tree)
    a2 <- treebalance::maxDelW(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::max_del_width(focal_tree),
                           treestats::max_del_width(ltab))
  }
})

test_that("normalization", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::max_del_width(focal_tree)
  c2 <- treestats::max_del_width(focal_tree, normalization = "tips")
  testthat::expect_lt(c2, c1)
  c3 <- treestats::max_del_width(treestats::phylo_to_l(focal_tree),
                            normalization = "tips")
  testthat::expect_equal(c2, c3)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::max_del_width(focal_tree))
    stats2 <- c(stats2, treestats::max_del_width(focal_tree,
                                                 normalization = "tips"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.5)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::max_del_width(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::max_del_width(list()),
    "input object has to be phylo or ltable"
  )
})
